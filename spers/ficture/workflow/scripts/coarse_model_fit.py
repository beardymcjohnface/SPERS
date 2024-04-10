import logging, pickle, sys
from joblib import Parallel, delayed
from sklearn.model_selection import train_test_split
from sklearn.utils import check_random_state
from scipy import sparse, stats
from sklearn.decomposition import LatentDirichletAllocation as LDA
import numpy as np
import pandas as pd


def read_input(infile):
    """
    Read in the transcripts TSV file into a pandas DF

    :param infile: filepath (str) of transcripts TSV
    :return: pandas dataframe, expect ["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]
    """
    if infile.endswith(".gz"):
        df = pd.read_csv(infile, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(infile, sep="\t")
    return df


def calculate_gene_weights(df):
    """
    Calculate weighting for all genes
    :param df: pandas dataframe long format
    :return: pandas dataframe with ["gene", "Count", "Weight"]
    """
    gene_weights = df.groupby(by=["gene"]).agg({"Count": sum}).reset_index()
    gene_weights.Weight = gene_weights.Count * 1. / gene_weights.Count.sum()
    return gene_weights


def df_to_mtx(df):
    """
    Convert pandas dataframe to matrix

    :param df: pandas input long dataframe ["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]
    :return: pandas wide dataframe [cols=gene, rows=hex_id, values=sum(Count)]
    """
    # Add hex bin IDs
    df["hex_id"] = df["xbin"].astype(str) + ":" + df["ybin"].astype(str)

    # Drop unneeded columns
    df = df[["hex_id", "gene", "Count"]]

    # Aggregate counts
    df = df.groupby(["hex_id", "gene"]).sum()

    # Pivot wide
    df = df.pivot_table(index="hex_id", columns="gene", values="Count", aggfunc="sum", fill_value=0)

    # Return dataframe
    return df


def gen_even_slices(n, n_packs):
    start = 0
    if n_packs < 1:
        raise ValueError("gen_even_slices got n_packs=%s, must be >=1" % n_packs)
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            yield np.arange(start, end)
            start = end


def chisq(k,info,total_k,total_umi):
    res = []
    if total_k <= 0:
        return res
    for name, v in info.iterrows():
        if v[k] <= 0:
            continue
        tab=np.zeros((2,2))
        tab[0,0]=v[str(k)]
        tab[0,1]=v["gene_total"]-tab[0,0]
        tab[1,0]=total_k-tab[0,0]
        tab[1,1]=total_umi-total_k-v["gene_total"]+tab[0,0]
        fd=tab[0,0]/total_k/tab[0,1]*(total_umi-total_k)
        if fd < 1:
            continue
        tab = np.around(tab, 0).astype(int) + 1
        chi2, p, dof, ex = stats.chi2_contingency(tab, correction=False)
        res.append([name,k,chi2,p,fd,v["gene_total"]])
    return res


def train_select_lda(mtx, hex_bins, **params):
    """
    Train and select lda

    :param mtx: pandas dataframe wide matrix (cols=gene, rows=hex_id)
    :param hex_bins: pandas daraframe long format from input
    :param params: parameters passed from config
    :return:
    """

    # convert matrix to CSR
    mtx_csr = sparse.coo_array(mtx).tocsr()

    # split into test and train
    train_mtx, test_mtx = train_test_split(mtx_csr, test_size=params["test_split"])
    n_train, _ = train_mtx.shape
    n_test, _ = test_mtx.shape

    # misc params
    factor_header = list(np.arange(params["n_components"]).astype(str))
    gene_index = {x:i for i,x in enumerate(list(mtx.columns))}
    coherence_scores = []
    model_results = {}

    # Iterate run LDA
    for r in range(params["iterations"]):
        # Initialise the model
        model = LDA(
            n_components=params["n_components"],
            learning_method=params["learning_method"],
            batch_size=params["batch_size"],
            total_samples=params["total_samples"],
            learning_offset=params["learning_offset"],
            learning_decay=params["learning_decay"],
            doc_topic_prior=params["doc_topic_prior"],
            n_jobs=params["threads"],
            verbose=0,
            random_state=params["random_state"])

        # Shuffle and score
        # params["rng"].shuffle(train_mtx) # todo: fix
        _ = model.partial_fit(train_mtx)
        score_train = model.score(train_mtx) / n_train
        score_test = model.score(test_mtx) / n_test

        # Report
        logging.debug(f"{r}: {score_train:.2f}, {score_test:.2f}")

        # Transform the test set
        test_mtx_transform = model.transform(test_mtx)

        # Get DE genes from the test data
        info = test_mtx.tocsc().T @ test_mtx_transform
        info = pd.DataFrame(info, columns=factor_header)

        info.index = list(mtx.columns)
        info["gene_total"] = info[factor_header].sum(axis=1)
        info.drop(index=info.index[info.gene_total < params["min_transcripts_scored"]], inplace=True)
        total_k = np.array(info[factor_header].sum(axis=0))
        total_umi = info[factor_header].sum().sum()
        res = []

        for k, kname in enumerate(factor_header):
            idx_slices = [idx for idx in gen_even_slices(len(info), params["threads"])]
            with Parallel(n_jobs=params["threads"], verbose=0) as parallel:
                result = parallel(
                    delayed(chisq)(
                        kname,info.iloc[idx, :].loc[:, [kname, "gene_total"]],total_k[k], total_umi
                    ) for idx in idx_slices
                )
            res += [item for sublist in result for item in sublist]

        chidf = pd.DataFrame(res, columns=["gene", "factor", "Chi2", "pval", "FoldChange", "gene_total"])
        chidf["Rank"] = chidf.groupby(by="factor")["Chi2"].rank(ascending=False)
        chidf.gene_total = chidf.gene_total.astype(int)
        chidf.sort_values(by=["factor", "Chi2"], ascending=[True, False], inplace=True)

        # Need the gene weights
        logging.debug("Calculating gene weights")
        gene_weights = calculate_gene_weights(hex_bins)

        # Compute a "coherence" score using top DE gene co-occurrence
        score = []
        for k in range(params["n_components"]):
            wd_idx = chidf.loc[chidf.factor.eq(str(k))].gene.iloc[:params["top_n_models"]].map(gene_index).values
            wd_idx = sorted(list(wd_idx), key=lambda x: -gene_weights.Weight.values[x])
            s = 0
            for ii in range(params["top_n_models"] - 1):
                for jj in range(ii + 1, params["top_n_models"]):
                    i = wd_idx[ii]
                    j = wd_idx[jj]
                    idx = mtx.indices[mtx.indptr[i]:mtx.indptr[i + 1]]
                    denom = mtx[:, [i]].toarray()[idx] * gene_weights.Weight.values[j] / gene_weights.Weight.values[i]
                    num = mtx[:, [j]].toarray()[idx]
                    s += (test_mtx_transform[idx, k].reshape((-1, 1)) * np.log(num / denom + 1)).sum()
            s0 = s / test_mtx_transform[:, k].sum()
            coherence_scores.append([r, k, s, s0])
            score.append(s0)

        model_results[r] = {"score_train": score_train, "score_test": score_test, "model": model, "coherence": score}

    # Save results
    pickle.dump(model_results, open(params["out_res"], "wb"))
    coherence_scores = pd.DataFrame(coherence_scores, columns=["R", "K", "Score0", "Score"])
    coherence_scores.to_csv(params["out_coh"], sep="\t", index=False)
    coherence_scores = coherence_scores.groupby(by="R").Score.mean()
    coherence_scores = coherence_scores.sort_values(ascending=False)
    best_model = model_results[coherence_scores.index[0]]["model"]

    # Minibatch index
    batch_index = hex_bins[["hex_id", "random_index"]]

    # refine with minibatches
    for minibatch_id in np.unique(batch_index.random_index):
        batch_hex_ids = batch_index[batch_index.random_index==minibatch_id].hex_id
        batch_mtx = sparse.coo_array(mtx[mtx.hex_id.isin(batch_hex_ids)]).tocsr()
        _ = best_model.partial_fit(batch_mtx)

    # Relabel factors
    weight = best_model.components_.sum(axis=1)
    ordered_k = np.argsort(weight)[::-1]
    best_model.components_ = best_model.components_[ordered_k, :]
    best_model.exp_dirichlet_component_ = best_model.exp_dirichlet_component_[ordered_k, :]

    # Rerun all units once and store results
    output_header = ["hex_id", "Count", "x", "y", "topK", "topP"] + factor_header
    dtp = {"topK": int, "Count": int, "unit": str}
    dtp.update({x: float for x in ["topP"] + factor_header})

    # Get final hex bin model scores
    mtx_csr_transform = best_model.transform(mtx_csr)

    # FIT OUTPUT
    fit_result = pd.DataFrame()
    fit_result["hex_id"] = mtx.index
    fit_result["Count"] = np.sum(mtx_csr, axis=1)

    # Merge hex bin coords
    hex_bin_coords = hex_bins[["hex_id","xbin","ybin"]].drop_duplicates().set_index("hex_id")
    fit_result = pd.concat([fit_result, hex_bin_coords], axis=1, join="inner")
    fit_result.rename(columns={"xbin":"x", "ybin":"y"})

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
    post_count["gene"] = mtx["gene"]
    post_count[["gene"] + factor_header].to_csv(
        params["out_pos"], sep="\t", index=False, float_format="%.2f", compression="gzip")

    # MODEL MATRIX OUTPUT
    best_model.feature_names_in_ = mtx["gene"]

    # best_model.log_norm_scaling_const_ = scale_const # todo update if log norm
    best_model.unit_sum_mean_ = np.mean(np.sum(mtx_csr, axis=1))

    # model matrix dataframe
    model_matrix = pd.DataFrame(best_model.components_.T, columns=factor_header, dtype="float64")
    model_matrix["gene"] = mtx["gene"]

    # write
    model_matrix[["gene"] + factor_header].to_csv(
        params["out_mtx"], sep="\t", index=False, float_format="%.4e", compression="gzip")

    # MODEL PICKLE OUTPUT
    pickle.dump(best_model, open(params["out_mdl"], "wb"))


def main(params=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)

    params["rng"] = check_random_state(params["random_state"])

    logging.debug("Reading in hex binned transcripts")
    hex_bins = read_input(kwargs["infile"])

    logging.debug("Reformatting the dataframe")
    hex_mtx = df_to_mtx(hex_bins)
    params["total_samples"], _ = hex_mtx.shape

    # TODO: add in options for log normalisation

    logging.debug("Iterative running LatentDirichletAllocation")
    train_select_lda(hex_mtx, hex_bins, **kwargs, **params)




if __name__ == "__main__":
    main(
        infile=snakemake.input[0],
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