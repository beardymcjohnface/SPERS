import logging, pickle, sys
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
from spers.ficture.workflow.scripts.generate_lda_model import df_to_mtx


def score_batch(df, lda_model, hex_width, xy_offsets):

    # Generate hex bins for offset combination
    logging.debug("Generating hex bins  (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_df = transcript_to_hex_bins(df, x_offset=xy_offsets[0], y_offset=xy_offsets[1], hex_width=hex_width)

    # Transform for scoring
    logging.debug("Transforming         (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_mtx = df_to_mtx(hex_df)
    hex_mtx = hex_mtx[lda_model.feature_names_in_]

    # Score bins
    logging.debug("Scoring              (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_transform = lda_model.transform(hex_mtx)

    logging.debug("Slicing best         (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    itr_result = pd.DataFrame({
        "hex_id": hex_mtx.index,
        "topK": np.argmax(hex_transform, axis=1),
        "topP": hex_transform.max(axis=1)
    })

    # Join transcript IDs
    logging.debug("Merging output       (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    hex_df = hex_df.merge(itr_result, on="hex_id", how="inner")

    logging.debug("Returning results    (x: " + str(xy_offsets[0]) + " and y: " + str(xy_offsets[1]) + ")")
    return hex_df[["transcript_id", "topK", "topP"]]


def iterate_and_score(df, lda_model, hex_width=None, offset_steps=None, step_size=None, threads=1, **params):
    """
    Generate offset hex bins and score with model

    :param df: pandas dataframe ["transcript_id", "x", "y", "gene", ...]
    :param lda_model: pickle of trained LDA model
    :param hex_width: int width of hex bins
    :param offset_steps: list int offset X/Y steps to perform
    :param step_size: int size of X/Y offset steps
    :param params:
    :return: pandas dataframe ["transcript_id", "best_hex_id", "best_topK", "best_topP"]
    """

    # Iterate offset steps and score
    all_results = Parallel(n_jobs=threads)(
        delayed(score_batch)(df, lda_model, hex_width, xy)
            for xy in ((x * step_size, y * step_size) for x in offset_steps for y in offset_steps)
    )

    # get top score for each transcript
    all_results = pd.concat(all_results).sort_values(
        by="topP", ascending=False).groupby("transcript_id").head(1).reset_index()

    # Done!
    return all_results


def main(in_tsv=None, in_mdl=None, out_tsv=None, log_file=None, threads=1, params=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    logging.debug("Running transcript_rescore.py")

    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_csv(in_tsv, sep="\t", compression="gzip")

    logging.debug("Reading in trained LDA model")
    lda_model = pickle.load(open(in_mdl, "rb"))
    lda_model.n_jobs = threads

    logging.debug("Generating overlapping hex bins and scoring")
    transcript_scores = iterate_and_score(transcripts_df, lda_model, threads=threads, **params)

    logging.debug("Writing best scores for each transcript")
    transcript_scores.to_csv(out_tsv, sep="\t", compression="gzip", index=False, float_format="%.2f")


if __name__ == "__main__":
    main(
        in_tsv=snakemake.input.tsv,
        in_mdl=snakemake.input.mdl,
        out_tsv=snakemake.output[0],
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        params=snakemake.params.params
    )