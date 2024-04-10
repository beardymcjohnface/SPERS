import logging, math
import numpy as np
import pandas as pd


def read_input(infile):
    """
    Read in the transcripts TSV file into a pandas DF

    :param infile: filepath (str) of transcripts TSV
    :return: pandas dataframe
    """
    if infile.endswith(".gz"):
        df = pd.read_csv(infile, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(infile, sep="\t")
    return df


def filter_min_transcripts_gene(df, min_transcripts_per_gene=None, min_transcripts_per_minibatch=None):
    """
    Filter transcripts for genes with total count < min gene count

    :param df: Pandas dataframe with ["random_index", "X", "Y", "gene", "Count"]
    :param min_transcripts_per_gene: int min required transcripts per gene
    :param min_transcripts_per_minibatch: int min required transcripts per minibatch - not implemented
    :return: Filtered Pandas dataframe with ["random_index", "X", "Y", "gene", "Count"]
    """

    # Collect total gene counts
    gene_counts = df.groupby(by=["gene"]).agg({"Count":sum}).reset_index()

    # List of genes above min cutoff
    gene_filter = list(gene_counts.loc[gene_counts["Count"]>=min_transcripts_per_gene]["gene"])

    # return the filtered dataframe
    return df.loc[df["gene"].isin(gene_filter)]


def transcript_to_hex_bins(df, scale=None):
    """
    Bin transcripts into hexagon bins based on param dimensions

    :param df: Pandas dataframe with ["random_index", "X", "Y", "gene", "Count"]
    :param scale: um size of hex bins
    :return: Pandas dataframe with ["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]
    """
    # calculate scale for Y
    scale_sqrt3 = scale * math.sqrt(3)

    # Scale coords
    df["xdiv"] = np.divmod(df["X"], scale / 2)[0]
    df["ydiv"] = np.divmod(df["Y"], scale_sqrt3 / 2)[0]

    # Nearest pair for X
    df["xn"] = scale / 2 * (df["xdiv"] + np.where(df["xdiv"] % 2 == 1, 1, 0))
    df["xns"] = scale / 2 * (df["xdiv"] + np.where(df["xdiv"] % 2 == 0, 1, 0))

    # Nearest pair for Y
    df["yn"] = scale_sqrt3 / 2 * (df["ydiv"] + np.where(df["ydiv"] % 2 == 1, 1, 0))
    df["yns"] = scale_sqrt3 / 2 * (df["ydiv"] + np.where(df["ydiv"] % 2 == 0, 1, 0))

    # Distances for each nearest pair
    df.loc[:, "d1"] = np.sqrt(
        (df["X"] - df["xn"]) * (df["X"] - df["xn"]) + (df["Y"] - df["yn"]) * (df["Y"] - df["yn"])
    )
    df.loc[:, "d2"] = np.sqrt(
        (df["X"] - df["xns"]) * (df["X"] - df["xns"]) + (df["Y"] - df["yns"]) * (df["Y"] - df["yns"])
    )

    # Retain the nearest coords
    df["xbin"] = np.where(df["d1"] > df["d2"], df["xn"], df["xns"])
    df["ybin"] = np.where(df["d1"] > df["d2"], df["yn"], df["yns"])

    # Drop unused columns and return
    return df[["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]]


def filter_bins_min_count(df, min_transcripts_per_hex=None):
    """
    Filter only hex bins with count > min_count

    :param df: pandas dataframe with ["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]
    :param min_transcripts_per_hex: int minimum "Count" transcripts per hex bin
    :return: pandas dataframe filtered with ["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]
    """
    # create hex bin IDs for counting and filtering
    df["hexid"] = df["xbin"].astype(str) + ":" + df["ybin"].astype(str)

    # Collect list of hex IDs with sum Count > min_count
    hex_filter = pd.DataFrame(df.groupby(["hexid"])["Count"].sum())
    hex_filter = list(hex_filter[hex_filter["Count"]>min_transcripts_per_hex].index)

    # Filter the dataframe
    df = df[df["hexid"].isin(hex_filter)]

    # Drop unused columns and return
    return df[["random_index", "X", "Y", "gene", "Count", "xbin", "ybin"]]


def main(infile=None, outfile=None, log_file=None, params=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading input transcripts: " + str(infile))
    transcripts_df = read_input(infile)

    logging.debug("Filtering low count genes")
    transcripts_df = filter_min_transcripts_gene(transcripts_df, min_transcripts_per_gene=params["min_transcripts_per_gene"])

    logging.debug("Hex binning the transcripts")
    transcripts_df = transcript_to_hex_bins(transcripts_df, scale=params["width"])

    logging.debug("Filtering hex bins")
    transcripts_df = filter_bins_min_count(transcripts_df, min_transcripts_per_hex=params["min_transcripts_per_hex"])

    logging.debug("Writing hex-binned transcripts")
    transcripts_df.to_csv(outfile, sep="\t", compression="gzip", index=False, float_format="%.2f")


if __name__ == "__main__":
    main(
        infile=snakemake.input[0],
        outfile=snakemake.output[0],
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )