import logging, sys
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins


def plot_transcripts(plot_df, out_png, **params):

    # Get scale for output plot
    x_scale = (plot_df["x"].max() - plot_df["x"].min()) / params["microns_per_inch"]
    y_scale = (plot_df["y"].max() - plot_df["y"].min()) / params["microns_per_inch"]
    plt.figure(figsize=(x_scale, y_scale))

    # Plot
    seaborn.set(
        rc={'axes.facecolor': params["background"], 'figure.facecolor': params["background"], 'axes.grid': False},
        font_scale=params["font_scale"])

    coarse_plot = seaborn.scatterplot(
        shuffle(plot_df).groupby("hex_id").head(params["transcripts_per_hex"]),
        x="x",
        y="y",
        hue="topK",
        s=params["point_scale"],
        palette=params["palette"],
        edgecolor=None)

    coarse_plot.legend([],[], frameon=False)
    coarse_plot.set_xlabel("X (microns)")
    coarse_plot.set_ylabel("Y (microns)")

    plt.savefig(out_png)
    # plt.savefig(out_svg)


def main(in_scr=None, in_trn=None, out_png=None, log_file=None, params=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_csv(in_trn, sep="\t", compression="gzip", usecols=["transcript_id", "x", "y"])

    logging.debug("Reading in transcript scores")
    scores_df = pd.read_csv(in_scr, sep="\t", compression="gzip", usecols=["transcript_id", "best_topK"])

    logging.debug("Joining transcript classifications")
    scores_df.columns = ["transcript_id", "topK"]
    transcripts_df = transcripts_df.merge(scores_df, on="transcript_id", how="inner")

    logging.debug("Calculating new hex bins for plotting")
    transcripts_df = transcript_to_hex_bins(transcripts_df, hex_width=params["hex_width"])

    logging.debug("Filtering low transcript bins")
    hex_counts = transcripts_df.groupby("hex_id").size()
    hex_filter = {k for k,v in hex_counts.items() if v >= params["min_transcripts_per_hex"]}
    transcripts_df = transcripts_df[transcripts_df["hex_id"].isin(hex_filter)]

    logging.debug("Plotting scored transcripts")
    plot_transcripts(transcripts_df, out_png, **params)


if __name__ == "__main__":
    main(
        in_scr=snakemake.input.scr,
        in_trn=snakemake.input.trn,
        out_png=snakemake.output.png,
        # out_svg=snakemake.output.svg,
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )