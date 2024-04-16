import logging
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
from sklearn.utils import shuffle


def plot_coarse_hex(fit_df, hex_df, transcripts_df, out_png, **params):
    """
    Generate a seaborn plot object for coarse scored transcripts
    :param fit_df: pandas dataframe of fit results ["hex_id", "topK"]
    :param hex_df: pandas dataframe of coarse hex bins ["hex_id", "transcript_id"]
    :param transcripts_df: pandas dataframe ["transcript_id", "x", "y"]
    :param params: config params passed from snakemake rule
    :return: plot object to save
    """

    # Merge the relevant columns
    plot_df = hex_df.merge(fit_df[["hex_id", "topK"]], on="hex_id", how="inner")
    plot_df = plot_df.merge(transcripts_df, on="transcript_id", how="inner")

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


def main(params=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)

    logging.debug("Reading in coarse model fit scores")
    fit_df = pd.read_csv(kwargs["in_fit"], sep="\t", compression="gzip", usecols=["hex_id", "topK", "x", "y"])

    logging.debug("Reading in coarse hex transcripts")
    hex_df = pd.read_csv(kwargs["in_hex"], sep="\t", compression="gzip", usecols=["hex_id", "transcript_id"])

    logging.debug("Reading in transcript coords")
    transcript_df = pd.read_csv(kwargs["in_trn"], sep="\t", compression="gzip", usecols=["transcript_id", "x", "y"])

    logging.debug("Plotting coarse-scored transcripts")
    plot_coarse_hex(fit_df, hex_df, transcript_df, kwargs["out_png"], **params)


if __name__ == "__main__":
    main(
        in_fit=snakemake.input.fit,
        in_hex=snakemake.input.hex,
        in_trn=snakemake.input.trn,
        out_png=snakemake.output.png,
        # out_svg=snakemake.output.svg,
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        params=snakemake.params.params
    )