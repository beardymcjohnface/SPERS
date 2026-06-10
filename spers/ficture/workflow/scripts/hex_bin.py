import logging, math
import numpy as np
import pandas as pd


def filter_min_transcripts_gene(df, min_transcripts_per_gene=None):
    """
    Filter transcripts for genes with total count < min gene count

    :param df: Pandas dataframe with ["transcript_id", "x", "y", "gene"]
    :param min_transcripts_per_gene: int min required transcripts per gene
    :return: Filtered Pandas dataframe with ["transcript_id", "x", "y", "gene"]
    """

    # Collect total gene counts (weight by "count" column if present, e.g. Visium HD)
    if "count" in df.columns:
        gene_counts = df.groupby("gene")["count"].sum()
    else:
        gene_counts = df.groupby("gene").size()

    # Vectorised membership: map each row's gene to its total and threshold
    return df[df["gene"].map(gene_counts) >= min_transcripts_per_gene]


def transcript_to_hex_bins(df, x_offset=0, y_offset=0, hex_width=None, **params):
    """
    Bin transcripts into hexagon bins based on param dimensions

    :param df: Pandas dataframe with ["x", "y", ...]
    :param scale: um size of hex bins
    :return: Pandas dataframe with ["hex_id", "xbin", "ybin", "x", "y", ...]
    """
    # float32 halves the memory of every intermediate array below; coordinate
    # precision in float32 (~1e-3 um at slide scale) is far finer than hex_width.
    half_w = np.float32(hex_width / 2)
    half_h = np.float32(hex_width * math.sqrt(3) / 2)

    # Coords as float32 with the offset folded in once (single add per axis)
    logging.debug("Init coords")
    xo = df["x"].to_numpy(dtype=np.float32) + np.float32(x_offset)
    yo = df["y"].to_numpy(dtype=np.float32) + np.float32(y_offset)

    # Lattice cell indices (integer-valued floats)
    logging.debug("Scale coords")
    x_div = np.floor_divide(xo, half_w)
    y_div = np.floor_divide(yo, half_h)

    # Parity is 0.0/1.0; one candidate centre adds it, the other its complement
    xpar = np.mod(x_div, 2)
    ypar = np.mod(y_div, 2)

    logging.debug("Nearest centre candidates")
    xn1 = half_w * (x_div + xpar)
    xn2 = half_w * (x_div + (1.0 - xpar))
    yn1 = half_h * (y_div + ypar)
    yn2 = half_h * (y_div + (1.0 - ypar))

    # Squared distances are sufficient for the nearest comparison (no sqrt),
    # and each delta is computed only once.
    logging.debug("Calculating distances")
    dx1 = xo - xn1; dy1 = yo - yn1
    dx2 = xo - xn2; dy2 = yo - yn2
    nearest1 = (dx1 * dx1 + dy1 * dy1) < (dx2 * dx2 + dy2 * dy2)
    del dx1, dy1, dx2, dy2

    # Retain the nearest centre coords
    logging.debug("Saving nearest centroids")
    xbin = np.where(nearest1, xn1, xn2)
    ybin = np.where(nearest1, yn1, yn2)
    df["xbin"] = xbin
    df["ybin"] = ybin

    # Build hex IDs from the chosen integer lattice indices. Factorising a single
    # integer key is much faster (and lower memory) than grouping two float columns.
    logging.debug("Creating hex IDs")
    xi = np.rint(xbin / half_w).astype(np.int64)
    yi = np.rint(ybin / half_h).astype(np.int64)
    yi0 = yi.min()
    key = (xi - xi.min()) * (yi.max() - yi0 + 1) + (yi - yi0)
    df["hex_id"] = pd.factorize(key, sort=False)[0]

    logging.debug("Done")
    return df


def filter_bins_min_count(df, min_transcripts_per_hex=None):
    """
    Filter only hex bins with count > min_count

    :param df: pandas dataframe with ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    :param min_transcripts_per_hex: int minimum "Count" transcripts per hex bin
    :return: pandas dataframe filtered with ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    """

    # Collect total hex transcript counts (weight by "count" column if present)
    if "count" in df.columns:
        hex_counts = df.groupby("hex_id")["count"].sum()
    else:
        hex_counts = df.groupby("hex_id").size()

    # Vectorised membership: map each row's hex to its total and threshold
    return df[df["hex_id"].map(hex_counts) >= min_transcripts_per_hex]
