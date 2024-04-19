import logging, pickle, sys
import pandas as pd
import numpy as np
from scipy import sparse

from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
from spers.ficture.workflow.scripts.coarse_model_fit import df_to_mtx


def iterate_and_score(df, lda_model, hex_width=None, offset_steps=None, step_size=None, **params):
    """
    Generate offset hex bins and score with model

    :param df: pandas dataframe ["transcript_id", "x", "y", "gene"]
    :param lda_model: pickle of trained LDA model
    :param hex_width: int width of hex bins
    :param offset_steps: list int offset X/Y steps to perform
    :param step_size: int size of X/Y offset steps
    :param params:
    :return: pandas dataframe ["transcript_id", "best_hex_id", "best_topK", "best_topP"]
    """

    # Initialise output
    output_columns = ["transcript_id", "best_topK", "best_topP", "best_hex_id"]
    out_df = df[["transcript_id"]]
    out_df["best_hex_id"] = "NA"
    out_df["best_topK"] = 0
    out_df["best_topP"] = 0

    # Iterate offset steps
    for x_step in offset_steps:
        x_offset = x_step * step_size
        for y_step in offset_steps:
            y_offset = y_step * step_size

            # Generate hex bins for offset combination
            hex_df = transcript_to_hex_bins(df, x_offset=x_offset, y_offset=y_offset, hex_width=hex_width)

            # Transform for scoring
            hex_mtx = df_to_mtx(hex_df)
            hex_csr = sparse.coo_array(hex_mtx[lda_model.feature_names_in_]).tocsr()

            # Score bins
            hex_transform = lda_model.transform(hex_csr)
            itr_result = pd.DataFrame({
                "hex_id": hex_mtx.index,
                "topK": np.argmax(hex_transform, axis=1),
                "topP": hex_transform.max(axis=1)
            })

            # Join transcript IDs
            hex_df = hex_df.merge(itr_result, on="hex_id", how="inner")

            # Join running best results
            out_df = out_df[output_columns].merge(hex_df[["hex_id", "transcript_id", "topK", "topP"]], on="transcript_id", how="outer")

            # Keep the running best
            out_df["best_topP"] = np.where(out_df["topK"] > out_df["best_topK"], out_df["topP"], out_df["best_topP"])
            out_df["best_hex_id"] = np.where(out_df["topK"] > out_df["best_topK"], out_df["hex_id"], out_df["best_hex_id"])
            out_df["best_topK"] = np.where(out_df["topK"] > out_df["best_topK"], out_df["topK"], out_df["best_topK"])

    # Done!
    return out_df[output_columns]


def main(in_tsv=None, in_mdl=None, out_tsv=None, log_file=None, threads=1, params=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    logging.debug("Running transcript_rescore.py")

    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_csv(in_tsv, sep="\t", compression="gzip")

    logging.debug("Reading in trained LDA model")
    lda_model = pickle.load(open(in_mdl, "rb"))
    lda_model.n_jobs = threads

    logging.debug("Generating overlapping hex bins and scoring")
    transcript_scores = iterate_and_score(transcripts_df, lda_model, **params)

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