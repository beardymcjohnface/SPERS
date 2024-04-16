import logging, math, sys
import numpy as np
import pandas as pd


def filter_min_transcripts_gene(df, min_transcripts_per_gene=None):
    """
    Filter transcripts for genes with total count < min gene count

    :param df: Pandas dataframe with ["transcript_id", "x", "y", "gene"]
    :param min_transcripts_per_gene: int min required transcripts per gene
    :return: Filtered Pandas dataframe with ["transcript_id", "x", "y", "gene"]
    """

    # Collect total gene counts
    gene_counts = df.groupby("gene").size()

    # List of genes above min cutoff
    gene_filter = {k for k,v in gene_counts.items() if v >= min_transcripts_per_gene}

    # return the filtered dataframe
    return df.loc[df["gene"].isin(gene_filter)]


def transcript_to_hex_bins(df, x_offset=0, y_offset=0, hex_width=None, **params):
    """
    Bin transcripts into hexagon bins based on param dimensions

    :param df: Pandas dataframe with ["x", "y", ...]
    :param scale: um size of hex bins
    :return: Pandas dataframe with ["hex_id", "xbin", "ybin", "x", "y", ...]
    """
    # calculate scale for y
    scale_sqrt3 = hex_width * math.sqrt(3)

    # numpy array
    logging.debug("Init numpy arrays")
    x_ar = np.array(df["x"])
    y_ar = np.array(df["y"])

    # Scale coords
    logging.debug("Scale coords")
    x_div = np.divmod(x_ar + x_offset, hex_width / 2)[0]
    y_div = np.divmod(y_ar + y_offset, scale_sqrt3 / 2)[0]

    # Nearest pair for x
    logging.debug("Nearest x coords")
    xn1 = hex_width / 2 * (x_div + np.where(x_div % 2 == 1, 1, 0))
    xn2 = hex_width / 2 * (x_div + np.where(x_div % 2 == 0, 1, 0))

    # Nearest pair for y
    logging.debug("Nearest y coords")
    yn1 = scale_sqrt3 / 2 * (y_div + np.where(y_div % 2 == 1, 1, 0))
    yn2 = scale_sqrt3 / 2 * (y_div + np.where(y_div % 2 == 0, 1, 0))

    # Distances for each nearest pair
    logging.debug("Calculating distances")
    dn1 = np.sqrt(
        (x_ar + x_offset - xn1) * (x_ar + x_offset - xn1) +
        (y_ar + y_offset - yn1) * (y_ar + y_offset - yn1)
    )
    dn2 = np.sqrt(
        (x_ar + x_offset - xn2) * (x_ar + x_offset - xn2) +
        (y_ar + y_offset - yn2) * (y_ar + y_offset - yn2)
    )

    # Retain the nearest coords
    logging.debug("Saving nearest centroids")
    xbin = np.where(dn1 < dn2, xn1, xn2)
    ybin = np.where(dn1 < dn2, yn1, yn2)

    # Return dataframe
    logging.debug("Adding to dataframe")
    df["xbin"] = xbin
    df["ybin"] = ybin

    # create hex bin IDs for counting and filtering
    logging.debug("Creating hex IDs")
    # TODO: abandon this format of hex_id? is slow
    df["hex_id"] = df["xbin"].astype(str) + "_" + df["ybin"].astype(str) + "_" + str(x_offset) + "_" + str(y_offset)

    logging.debug("Done")
    return df


def filter_bins_min_count(df, min_transcripts_per_hex=None):
    """
    Filter only hex bins with count > min_count

    :param df: pandas dataframe with ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    :param min_transcripts_per_hex: int minimum "Count" transcripts per hex bin
    :return: pandas dataframe filtered with ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    """

    # Collect total hex transcript counts
    hex_counts = df.groupby("hex_id").size()

    # List of hex_ids above min cutoff
    hex_filter = {k for k,v in hex_counts.items() if v >= min_transcripts_per_hex}

    # Filter the dataframe
    df = df[df["hex_id"].isin(hex_filter)]

    # Drop unused columns and return
    return df


def main(infile=None, outfile=None, log_file=None, params=None):
    file_handler = logging.FileHandler(filename=log_file)
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    handlers = [file_handler, stdout_handler]

    logging.basicConfig(handlers=handlers,level=logging.DEBUG)

    logging.debug("Reading input transcripts: " + str(infile))
    transcripts_df = pd.read_csv(infile, sep="\t", compression="gzip")

    if "min_transcripts_per_gene" in params.keys():
        logging.debug("Filtering low count genes")
        transcripts_df = filter_min_transcripts_gene(transcripts_df, min_transcripts_per_gene=params["min_transcripts_per_gene"])

    logging.debug("Hex binning the transcripts")
    transcripts_df = transcript_to_hex_bins(transcripts_df, **params)

    if "min_transcripts_per_hex" in params.keys():
        logging.debug("Filtering hex bins")
        transcripts_df = filter_bins_min_count(
            transcripts_df,
            min_transcripts_per_hex=params["min_transcripts_per_hex"]
        )

    logging.debug("Writing output file")
    # TODO: don't write the hex bins, just calculate them as needed, writing is the bottleneck
    # pickle.dump(transcripts_df, open(outfile, "wb"))
    transcripts_df.to_csv(outfile, sep="\t", index=False, float_format="%.2f", compression="gzip")


if __name__ == "__main__":
    main(
        infile=snakemake.input[0],
        outfile=snakemake.output[0],
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )