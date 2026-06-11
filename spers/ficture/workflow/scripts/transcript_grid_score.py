import logging, pickle, sys, warnings
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

# The model carries feature_names_in_ (used to align columns), but we score
# name-less sparse matrices, which triggers a benign sklearn warning per batch.
warnings.filterwarnings("ignore", message="X does not have valid feature names")

try:
    from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
    from spers.ficture.workflow.scripts.generate_lda_model import df_to_mtx, align_to_features
except ModuleNotFoundError:
    # Fall back to sibling imports when run as a Snakemake script in an
    # isolated conda env where the spers package is not installed.
    from hex_bin import transcript_to_hex_bins
    from generate_lda_model import df_to_mtx, align_to_features


def score_batch(df, lda_model, hex_width, xy_offsets):

    # Generate hex bins for offset combination
    logging.debug("Generating hex bins  (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_df = transcript_to_hex_bins(df, x_offset=xy_offsets[0], y_offset=xy_offsets[1], hex_width=hex_width)

    # Transform for scoring (sparse matrix, columns aligned to the model genes)
    logging.debug("Transforming         (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_mtx, hex_ids, hex_genes = df_to_mtx(hex_df)
    hex_mtx = align_to_features(hex_mtx, hex_genes, lda_model.feature_names_in_)

    # Score bins
    logging.debug("Scoring              (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_transform = lda_model.transform(hex_mtx)

    logging.debug("Slicing best         (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    # Keep the full per-factor score vector alongside the top factor/probability
    factor_header = [str(i) for i in range(hex_transform.shape[1])]
    itr_result = pd.DataFrame(hex_transform, columns=factor_header, index=pd.Index(hex_ids, name="hex_id"))
    itr_result["topK"] = np.argmax(hex_transform, axis=1)
    itr_result["topP"] = hex_transform.max(axis=1)
    itr_result = itr_result.reset_index()  # brings "hex_id" back as a column

    # Join transcript IDs
    logging.debug("Merging output       (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_df = hex_df.merge(itr_result, on="hex_id", how="inner")

    logging.debug("Returning results    (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    return hex_df[["transcript_id", "topK", "topP"] + factor_header]


def iterate_and_score(df, lda_model, hex_width=None, offset_steps=None, step_size=None, threads=1, **params):
    """
    Generate offset hex bins and score with model

    :param df: pandas dataframe ["transcript_id", "x", "y", "gene", ...]
    :param lda_model: pickle of trained LDA model
    :param hex_width: int width of hex bins
    :param offset_steps: list int offset X/Y steps to perform
    :param step_size: int size of X/Y offset steps
    :param params:
    :return: pandas dataframe ["transcript_id", "topK", "topP", "0", "1", ... per-factor scores]
    """

    # Iterate offset steps and score
    all_results = Parallel(n_jobs=threads)(
        delayed(score_batch)(df, lda_model, hex_width, xy)
            for xy in ((x * step_size, y * step_size) for x in offset_steps for y in offset_steps)
    )

    # get top score for each transcript (keeps that overlap's full factor vector)
    all_results = pd.concat(all_results).sort_values(
        by="topP", ascending=False).groupby("transcript_id").head(1).reset_index(drop=True)

    # Done!
    return all_results


def score_spot_level(df, lda_model, hex_width=None, **params):
    """
    Score spot-level hex bins directly, without expanding back to transcripts.

    For high-density platforms (e.g. Visium HD) the per-transcript output is huge
    and redundant (every transcript in a spot shares the same hex). This bins the
    transcripts once at the grid_score hex_width and returns one row per hex bin
    with its full factor-score vector.

    :param df: pandas dataframe ["transcript_id", "x", "y", "gene", "count", ...]
    :param lda_model: trained LDA model
    :param hex_width: int width of hex bins (= grid_score hex_width)
    :return: pandas dataframe ["hex_id", "x", "y", "Count", "topK", "topP", "0", "1", ...]
    """
    logging.debug("Binning transcripts into spot-level hex bins")
    hex_df = transcript_to_hex_bins(df, hex_width=hex_width)

    logging.debug("Building hex count matrix")
    hex_mtx, hex_ids, hex_genes = df_to_mtx(hex_df)
    hex_mtx = align_to_features(hex_mtx, hex_genes, lda_model.feature_names_in_)

    logging.debug("Scoring hex bins")
    hex_transform = lda_model.transform(hex_mtx)

    factor_header = [str(i) for i in range(hex_transform.shape[1])]
    result = pd.DataFrame(hex_transform, columns=factor_header, index=pd.Index(hex_ids, name="hex_id"))
    result["topK"] = np.argmax(hex_transform, axis=1)
    result["topP"] = hex_transform.max(axis=1)
    result["Count"] = np.asarray(hex_mtx.sum(axis=1)).ravel()

    # Attach hex-bin centroid coordinates
    coords = hex_df[["hex_id", "xbin", "ybin"]].drop_duplicates("hex_id").set_index("hex_id")
    result = result.join(coords).rename(columns={"xbin": "x", "ybin": "y"})
    result = result.reset_index()  # brings "hex_id" back as a column

    return result[["hex_id", "x", "y", "Count", "topK", "topP"] + factor_header]


def main(in_tsv=None, in_mdl=None, out_tsv=None, log_file=None, threads=1, params=None, platform=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    logging.debug("Running transcript_grid_score.py")

    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_pickle(in_tsv)

    logging.debug("Reading in trained LDA model")
    lda_model = pickle.load(open(in_mdl, "rb"))
    lda_model.n_jobs = threads

    if platform == "visiumhd":
        # Spot-level (hex-bin) scores: one row per hex bin, no transcript expansion
        logging.debug("Spot-level hex-bin scoring (visiumhd)")
        scores = score_spot_level(transcripts_df, lda_model, hex_width=params["hex_width"])
    else:
        logging.debug("Generating overlapping hex bins and scoring per transcript")
        scores = iterate_and_score(transcripts_df, lda_model, threads=threads, **params)

    logging.debug("Writing scores")
    scores.to_pickle(out_tsv)


if __name__ == "__main__":
    main(
        in_tsv=snakemake.input.tsv,
        in_mdl=snakemake.input.mdl,
        out_tsv=snakemake.output[0],
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        params=snakemake.params.params,
        platform=snakemake.params.platform
    )