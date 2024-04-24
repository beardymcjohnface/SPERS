import logging
import pandas as pd
from sklearn.utils import shuffle

from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
from spers.ficture.workflow.scripts.plot import plot_ficture


def main(in_scr=None, in_trn=None, out_png=None, log_file=None, params=None, plot=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_csv(in_trn, sep="\t", compression="gzip", usecols=["transcript_id", "x", "y"])

    logging.debug("Reading in transcript scores")
    scores_df = pd.read_csv(in_scr, sep="\t", compression="gzip", usecols=["transcript_id", "topK"])

    logging.debug("Joining transcript classifications")
    scores_df.columns = ["transcript_id", "topK"]
    transcripts_df = transcripts_df.merge(scores_df, on="transcript_id", how="inner")

    logging.debug("Calculating new hex bins for plotting")
    transcripts_df = transcript_to_hex_bins(transcripts_df, hex_width=params["grid_score"]["hex_width"])

    logging.debug("Filtering low transcript bins")
    hex_counts = transcripts_df.groupby("hex_id").size()
    hex_filter = {k for k,v in hex_counts.items() if v >= params["filter"]["min_transcripts_per_hex"]}
    transcripts_df = transcripts_df[transcripts_df["hex_id"].isin(hex_filter)]

    logging.debug("Downsampling transcripts for plotting")
    transcripts_df = shuffle(transcripts_df).groupby("hex_id").head(params["filter"]["max_transcripts_per_hex"])

    logging.debug("Plotting scored transcripts")
    plot_ficture(
        transcripts_df,
        out_png,
        point_scale = 0.00001,
        font_scale = 0.1,
        **plot)


if __name__ == "__main__":
    main(
        in_scr=snakemake.input.scr,
        in_trn=snakemake.input.trn,
        out_png=snakemake.output.png,
        log_file=snakemake.log[0],
        params=snakemake.params.params,
        plot=snakemake.params.plot
    )