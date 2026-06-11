import logging, pickle, sys
from joblib import Parallel, delayed, parallel_backend
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
from scipy import sparse, stats
from sklearn.decomposition import LatentDirichletAllocation as LDA
import numpy as np
import pandas as pd

try:
    # Installed-package import (interactive use / README examples)
    import spers.ficture.workflow.scripts.hex_bin as hex_bin
    import spers.ficture.workflow.scripts.spatial_minibatch as minibatch
except ModuleNotFoundError:
    # Snakemake adds the script's directory to sys.path, so the sibling
    # modules are importable even when the spers package is not installed
    # (e.g. inside an isolated conda env).
    import hex_bin
    import spatial_minibatch as minibatch


def calculate_gene_weights(df):
    """
    Calculate weighting for all genes
    :param df: pandas dataframe ["gene", ...]
    :return: pandas dataframe with ["gene", "Weight"]
    """
    if "count" in df.columns:
        gene_weights = df.groupby("gene")["count"].sum().reset_index()
    else:
        gene_weights = df.groupby("gene").size().reset_index()
    gene_weights.columns = ["gene", "Count"]
    gene_weights["Weight"] = gene_weights.Count * 1. / gene_weights.Count.sum()
    return gene_weights[["gene","Weight"]]


def df_to_mtx(df):
    """
    Build a sparse hex x gene count matrix directly (no dense pivot_table).

    Factorises hex_id and gene to integer row/column indices and constructs a
    scipy CSR matrix; duplicate (hex, gene) entries are summed on conversion, so
    no groupby is needed. This avoids materialising a dense hex x gene array
    (catastrophic at Visium HD gene counts) and the dense->sparse round-trip.

    :param df: pandas input long dataframe ["hex_id", "gene", ("count")]
    :return: (scipy.sparse.csr_matrix, hex_ids ndarray, genes ndarray)
             rows aligned to hex_ids, columns aligned to genes (both sorted)
    """
    rows, hex_ids = pd.factorize(df["hex_id"], sort=True)
    cols, genes = pd.factorize(df["gene"], sort=True)

    # Weight by "count" column if present (e.g. Visium HD), else one per row
    if "count" in df.columns:
        data = df["count"].to_numpy()
    else:
        data = np.ones(len(df), dtype=np.int64)

    mtx = sparse.coo_matrix((data, (rows, cols)), shape=(len(hex_ids), len(genes)))
    return mtx.tocsr(), np.asarray(hex_ids), np.asarray(genes)


def align_to_features(mtx_csr, genes, feature_names):
    """
    Reorder/subset sparse-matrix columns to match a model's feature order.

    Genes in feature_names that are absent from `genes` become all-zero columns,
    matching the previous DataFrame `mtx[feature_names]` behaviour but on a
    sparse matrix.

    :param mtx_csr: scipy CSR matrix (rows=hex, cols=genes)
    :param genes: sequence of gene labels for the matrix columns
    :param feature_names: target gene order (e.g. lda_model.feature_names_in_)
    :return: scipy CSR matrix with columns ordered as feature_names
    """
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    zero_col = len(genes)
    src = np.array([gene_to_idx.get(g, zero_col) for g in feature_names])

    # Append one zero column so missing genes map to it, then reindex columns
    padded = sparse.hstack(
        [mtx_csr, sparse.csr_matrix((mtx_csr.shape[0], 1), dtype=mtx_csr.dtype)],
        format="csc")
    return padded[:, src].tocsr()


def _vectorized_chisq(info, factor_header):
    """
    Vectorised 2x2 chi-square enrichment of each gene within each factor.

    Replaces the per-gene scipy.chi2_contingency loop with the closed-form 2x2
    statistic computed over the whole (gene x factor) matrix at once, which is
    far faster and removes the per-factor joblib pools. Returns the same columns
    as before: ["gene", "factor", "Chi2", "pval", "FoldChange", "gene_total", "Rank"].

    :param info: pandas dataframe indexed by gene, columns = factor_header (+ gene_total)
    :param factor_header: list of factor column names (strings)
    :return: pandas dataframe of per (gene, factor) enrichment stats, sorted
    """
    genes = list(info.index)
    a = info[factor_header].to_numpy(dtype="float64")        # gene-in-factor counts (G x K)
    gene_total = a.sum(axis=1)                               # (G,)
    total_k = a.sum(axis=0)                                  # (K,)
    total_umi = a.sum()

    b = gene_total[:, None] - a
    c = total_k[None, :] - a
    d = total_umi - total_k[None, :] - gene_total[:, None] + a

    # Enrichment fold change (pre-pseudocount), used to keep only enriched genes
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = (a / total_k[None, :]) / b * (total_umi - total_k[None, :])

    keep = (a > 0) & (fold >= 1)

    # Pseudocounted integer table (matches np.around(tab, 0).astype(int) + 1)
    ai = np.rint(a) + 1; bi = np.rint(b) + 1
    ci = np.rint(c) + 1; di = np.rint(d) + 1
    n = ai + bi + ci + di
    chi2 = n * (ai * di - bi * ci) ** 2 / ((ai + bi) * (ci + di) * (ai + ci) * (bi + di))
    pval = stats.chi2.sf(chi2, 1)

    gi, ki = np.nonzero(keep)
    chidf = pd.DataFrame({
        "gene": [genes[g] for g in gi],
        "factor": [factor_header[k] for k in ki],
        "Chi2": chi2[gi, ki],
        "pval": pval[gi, ki],
        "FoldChange": fold[gi, ki],
        "gene_total": gene_total[gi].astype(int),
    })
    chidf["Rank"] = chidf.groupby(by="factor")["Chi2"].rank(ascending=False)
    chidf.sort_values(by=["factor", "Chi2"], ascending=[True, False], inplace=True)
    return chidf


def _fit_one_model(r, train_mtx, test_mtx, test_mtx_csc, gene_weights, gene_index,
                   gene_names, factor_header, n_train, n_test, inner_jobs, params):
    """
    Fit a single candidate LDA model and score its coherence.

    Designed to run inside a joblib worker. Any spare cores (threads not consumed
    by the outer model pool) are given to the LDA via a threading backend so we
    don't nest process pools. Each model gets its own seed for independent fits.
    """
    seed = params["random_state"] + r
    model = LDA(**params["lda"], n_jobs=inner_jobs, verbose=0, random_state=seed)

    # Independent shuffle per model (own seed). Threading backend lets the LDA
    # E-step use the leftover cores without spawning nested processes.
    train_shuffled = shuffle(train_mtx, random_state=seed)
    with parallel_backend("threading", n_jobs=inner_jobs):
        model.partial_fit(train_shuffled)

        # model.score is two extra full passes and is not used for selection
        # (coherence is), so it is skipped unless explicitly requested.
        if params["train"].get("score_models", False):
            score_train = model.score(train_mtx) / n_train
            score_test = model.score(test_mtx) / n_test
        else:
            score_train = score_test = float("nan")

        test_mtx_transform = model.transform(test_mtx)
    logging.debug(f"{r}: {score_train:.2f}, {score_test:.2f}")

    # DE genes from the test data (vectorised chi-square)
    info = pd.DataFrame(test_mtx_csc.T @ test_mtx_transform, columns=factor_header, index=gene_names)
    info["gene_total"] = info[factor_header].sum(axis=1)
    info = info[info["gene_total"] >= params["train"]["min_transcripts_scored"]]
    chidf = _vectorized_chisq(info, factor_header)

    # Coherence score using top DE gene co-occurrence
    n_top = params["train"]["output_models"]
    weights = gene_weights.Weight.values
    score = []
    coherence_rows = []
    for k in range(params["lda"]["n_components"]):
        wd_idx = chidf.loc[chidf.factor.eq(str(k))].gene.iloc[:n_top].map(gene_index).values
        wd_idx = sorted(list(wd_idx), key=lambda x: -weights[x])
        s = 0
        for ii in range(n_top - 1):
            for jj in range(ii + 1, n_top):
                i = wd_idx[ii]
                j = wd_idx[jj]
                idx = test_mtx_csc.indices[test_mtx_csc.indptr[i]:test_mtx_csc.indptr[i + 1]]
                denom = test_mtx_csc[:, [i]].toarray()[idx] * weights[j] / weights[i]
                num = test_mtx_csc[:, [j]].toarray()[idx]
                s += (test_mtx_transform[idx, k].reshape((-1, 1)) * np.log(num / denom + 1)).sum()
        s0 = s / test_mtx_transform[:, k].sum()
        coherence_rows.append([r, k, s, s0])
        score.append(s0)

    return r, {"score_train": score_train, "score_test": score_test,
               "model": model, "coherence": score}, coherence_rows


def train_select_lda(mtx_csr, hex_ids, genes, transcripts_df, **params):
    """
    Train and select lda

    :param mtx_csr: scipy CSR hex x gene count matrix
    :param hex_ids: ndarray of hex_id row labels (aligned to mtx_csr rows)
    :param genes: ndarray of gene column labels (aligned to mtx_csr columns)
    :param transcripts_df: pandas long dataframe ["hex_id", "transcript_id", "xbin", "ybin", "gene", ...]
    :param params: parameters passed from config
    :return:
    """

    # split into test and train
    train_mtx, test_mtx = train_test_split(mtx_csr, test_size=params["train"]["test_split"])
    test_mtx_csc = test_mtx.tocsc()
    n_train, _ = train_mtx.shape
    n_test, _ = test_mtx.shape

    # Need the gene weights
    logging.debug("Calculating gene weights")
    gene_weights = calculate_gene_weights(transcripts_df)

    # misc params
    factor_header = list(np.arange(params["lda"]["n_components"]).astype(str))
    gene_names = list(genes)
    gene_index = {x: i for i, x in enumerate(gene_names)}

    # Train the independent candidate models in parallel. Cores left over after
    # the model pool (threads // n_workers) are handed to each LDA's E-step via a
    # threading backend, so spare CPUs are not left idle.
    n_workers = min(params["threads"], params["train"]["generate_models"])
    inner_jobs = max(1, params["threads"] // n_workers)
    logging.debug(
        f"Training {params['train']['generate_models']} models on {n_workers} worker(s) x {inner_jobs} thread(s)")
    fitted = Parallel(n_jobs=n_workers)(
        delayed(_fit_one_model)(
            r, train_mtx, test_mtx, test_mtx_csc, gene_weights, gene_index,
            gene_names, factor_header, n_train, n_test, inner_jobs, params)
        for r in range(params["train"]["generate_models"]))

    # Collect results
    coherence_scores = []
    model_results = {}
    for r, res, coh_rows in fitted:
        model_results[r] = res
        coherence_scores.extend(coh_rows)

    # Save results
    logging.debug("Saving model results and coherence scores")
    pickle.dump(model_results, open(params["out_res"], "wb"))
    coherence_scores = pd.DataFrame(coherence_scores, columns=["R", "K", "Score0", "Score"])
    coherence_scores.to_csv(params["out_coh"], sep="\t", index=False)
    coherence_scores = coherence_scores.groupby(by="R").Score.mean()
    coherence_scores = coherence_scores.sort_values(ascending=False)
    best_model = model_results[coherence_scores.index[0]]["model"]

    # refine with minibatches (integer row slicing on the sparse matrix instead
    # of an O(n_batches x n_hexes) isin scan + sparse rebuild per batch)
    logging.debug("Refining best model with minibatches")
    hexid_to_row = {h: i for i, h in enumerate(hex_ids)}
    for minibatch_hex_ids in minibatch.minibatch_transcripts(transcripts_df, **params["bin"]):
        rows = [hexid_to_row[h] for h in minibatch_hex_ids if h in hexid_to_row]
        if len(rows) > 1:
            _ = best_model.partial_fit(mtx_csr[rows])

    # Relabel factors
    weight = best_model.components_.sum(axis=1)
    ordered_k = np.argsort(weight)[::-1]
    best_model.components_ = best_model.components_[ordered_k, :]
    best_model.exp_dirichlet_component_ = best_model.exp_dirichlet_component_[ordered_k, :]

    # Rerun all units once and store results
    output_header = ["hex_id", "Count", "x", "y", "topK", "topP"] + factor_header
    dtp = {"topK": int, "Count": int, "hex_id": str}
    dtp.update({x: float for x in ["topP"] + factor_header})

    # Get final hex bin model scores
    mtx_csr_transform = best_model.transform(mtx_csr)

    # FIT OUTPUT
    fit_result = pd.DataFrame(index=pd.Index(hex_ids, name="hex_id"))
    fit_result["Count"] = np.asarray(mtx_csr.sum(axis=1)).ravel()

    # Merge hex bin coords
    hex_bin_coords = transcripts_df[["hex_id","xbin","ybin"]]
    hex_bin_coords = hex_bin_coords.drop_duplicates().set_index("hex_id")
    fit_result = pd.concat([fit_result, hex_bin_coords], axis=1, join="inner").reset_index()
    fit_result = fit_result.rename(columns={"xbin":"x", "ybin":"y"})

    # Top K/P
    fit_result["topK"] = np.argmax(mtx_csr_transform, axis=1).astype(int)
    fit_result["topP"] = np.max(mtx_csr_transform, axis=1)

    # Merge with model scores
    fit_result = pd.concat((fit_result, pd.DataFrame(mtx_csr_transform, columns=factor_header)), axis=1).astype(dtp)

    # write model fit
    fit_result[output_header].to_csv(params["out_fit"], sep="\t", float_format="%.4e", index=False, header=True, compression="gzip")

    # POSTERIOR COUNT OUTPUT
    post_count = np.array(mtx_csr_transform.T @ mtx_csr)
    post_count = pd.DataFrame(post_count.T, columns=factor_header, dtype="float64")
    post_count["gene"] = gene_weights["gene"]
    post_count[["gene"] + factor_header].to_csv(
        params["out_pos"], sep="\t", index=False, float_format="%.2f", compression="gzip")

    # MODEL MATRIX OUTPUT
    best_model.feature_names_in_ = gene_weights["gene"]

    # best_model.log_norm_scaling_const_ = scale_const # todo update if log norm
    best_model.unit_sum_mean_ = np.mean(np.asarray(mtx_csr.sum(axis=1)))

    # model matrix dataframe
    model_matrix = pd.DataFrame(best_model.components_.T, columns=factor_header, dtype="float64")
    model_matrix["gene"] = gene_weights["gene"]

    # write
    model_matrix[["gene"] + factor_header].to_csv(
        params["out_mtx"], sep="\t", index=False, float_format="%.4e", compression="gzip")

    # MODEL PICKLE OUTPUT
    pickle.dump(best_model, open(params["out_mdl"], "wb"))


def main(params=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    logging.debug("Reading transcripts")
    transcripts_df = pd.read_pickle(kwargs["in_tsv"])

    logging.debug("Filtering low count genes")
    transcripts_df = hex_bin.filter_min_transcripts_gene(transcripts_df, min_transcripts_per_gene=params["bin"]["min_transcripts_per_gene"])

    logging.debug("Calculating hex bins")
    transcripts_df = hex_bin.transcript_to_hex_bins(transcripts_df, x_offset=0, y_offset=0, hex_width=params["bin"]["hex_width"])

    logging.debug("Filtering low count hex bins")
    transcripts_df = hex_bin.filter_bins_min_count(transcripts_df, min_transcripts_per_hex=params["bin"]["min_transcripts_per_hex"])

    logging.debug("Initialising minibatch parameters")
    params["bin"] = minibatch.batch_dimensions(transcripts_df, **params["bin"])

    logging.debug("Reformatting the dataframe")
    mtx_csr, hex_ids, genes = df_to_mtx(transcripts_df)

    # TODO: add in options for log normalisation

    logging.debug("Iterative running LatentDirichletAllocation")
    train_select_lda(mtx_csr, hex_ids, genes, transcripts_df, **kwargs, **params)




if __name__ == "__main__":
    main(
        in_tsv=snakemake.input.tsv,
        out_fit=snakemake.output.fit,
        out_res=snakemake.output.res,
        out_coh=snakemake.output.coh,
        out_pos=snakemake.output.pos,
        out_mtx=snakemake.output.mtx,
        out_mdl=snakemake.output.mdl,
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        params=snakemake.params.params
    )