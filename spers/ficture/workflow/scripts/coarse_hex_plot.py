import logging
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
from sklearn.utils import shuffle


def plot_coarse_hex(fit_df, hex_df, out_png, **params):
    """
    Generate a seaborn plot object for coarse scored transcripts
    :param fit_df: pandas dataframe of fit results ["hex_id", "topK"]
    :param hex_df: pandas dataframe of coarse hex bins ["hex_id", "X", "Y", "gene", "Count", "xbin", "ybin"]
    :param params: config params passed from snakemake rule
    :return: plot object to save
    """
    hex_df = hex_df.merge(fit_df, on="hex_id", how="inner")

    x_span = hex_df["X"].max() - hex_df["X"].min()
    y_span = hex_df["Y"].max() - hex_df["Y"].min()

    if x_span > y_span:
        plt.figure(figsize=(params["size_scale"], params["size_scale"] * y_span / x_span))
    else:
        plt.figure(figsize=(params["size_scale"] * x_span / y_span, params["size_scale"]))

    seaborn.scatterplot(
        shuffle(hex_df).groupby("hex_id").head(params["transcripts_per_hex"]),
        x="X",
        y="Y",
        hue="topK",
        s=params["point_scale"],
        palette=params["palette"])

    plt.savefig(out_png)


def main(params=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)

    logging.debug("Reading in coarse model fit scores")
    fit_df = pd.read_csv(kwargs["in_fit"], sep="\t", compression="gzip")

    logging.debug("Reading in coarse hex transcripts")
    hex_df = pd.read_csv(kwargs["in_hex"], sep="\t", compression="gzip")

    logging.debug("Plotting coarse-scored transcripts")
    plot_coarse_hex(fit_df, hex_df, kwargs["out_plt"], **params)


if __name__ == "__main__":
    main(
        in_fit=snakemake.input.fit,
        in_hex=snakemake.input.hex,
        out_plt=snakemake.output[0],
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        params=snakemake.params.params
    )