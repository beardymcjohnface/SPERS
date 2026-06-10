import gzip
import json
import logging

import numpy as np
import pandas as pd


def read_microns_per_pixel(scalefactors_json):
    """
    Read the micron-to-pixel scale factor from a Space Ranger scalefactors json.

    :param scalefactors_json: str (path) to scalefactors_json.json
    :return: float microns per full-resolution pixel
    """
    with open(scalefactors_json) as f:
        scalefactors = json.load(f)

    for key in ("microns_per_pixel", "microns_per_pixel_002um"):
        if key in scalefactors:
            return float(scalefactors[key])

    raise KeyError(
        "Could not find 'microns_per_pixel' in " + str(scalefactors_json)
    )


def read_barcodes(barcodes_tsv):
    """
    Read the ordered barcode list (matrix column order, 1-based in the mtx).

    :param barcodes_tsv: str (path) barcodes.tsv.gz
    :return: numpy array of barcode strings (0-based index)
    """
    return pd.read_csv(
        barcodes_tsv, sep="\t", header=None, compression="gzip"
    )[0].to_numpy()


def read_features(features_tsv):
    """
    Read the ordered feature list (matrix row order, 1-based in the mtx).

    Columns are: gene_id, gene (symbol), feature_type.

    :param features_tsv: str (path) features.tsv.gz
    :return: pandas dataframe ["gene_id", "gene", "feature_type"]
    """
    features = pd.read_csv(features_tsv, sep="\t", header=None, compression="gzip")
    names = ["gene_id", "gene", "feature_type"]
    features.columns = names[: features.shape[1]]
    return features


def read_positions(positions_parquet):
    """
    Read the per-barcode spatial coordinates.

    :param positions_parquet: str (path) tissue_positions.parquet
    :return: pandas dataframe indexed by "barcode" with pixel row/col columns
    """
    positions = pd.read_parquet(positions_parquet)
    # Space Ranger uses "barcode" as the key column
    positions = positions.set_index("barcode")
    return positions[["pxl_row_in_fullres", "pxl_col_in_fullres"]]


def count_mtx_header_lines(matrix_mtx):
    """
    Count the leading MatrixMarket comment ('%') lines so we can skip them
    plus the single size line that follows.

    :param matrix_mtx: str (path) matrix.mtx.gz
    :return: (int n_comment_lines, tuple size_line (n_features, n_barcodes, nnz))
    """
    n_comments = 0
    with gzip.open(matrix_mtx, "rt") as f:
        for line in f:
            if line.startswith("%"):
                n_comments += 1
            else:
                size = tuple(int(x) for x in line.split())
                return n_comments, size
    raise ValueError("No size line found in " + str(matrix_mtx))


def read_matrix(matrix_mtx, n_skip):
    """
    Read the MatrixMarket coordinate body as a long dataframe.

    :param matrix_mtx: str (path) matrix.mtx.gz
    :param n_skip: int number of header lines to skip (comments + size line)
    :return: pandas dataframe ["feature_idx", "barcode_idx", "count"] (1-based idx)
    """
    return pd.read_csv(
        matrix_mtx,
        compression="gzip",
        sep=r"\s+",
        skiprows=n_skip,
        header=None,
        names=["feature_idx", "barcode_idx", "count"],
        dtype={"feature_idx": np.int32, "barcode_idx": np.int32, "count": np.int32},
    )


def visiumhd_to_transcripts(
    matrix_mtx, barcodes_tsv, features_tsv, positions_parquet, scalefactors_json, params
):
    """
    Convert a Visium HD binned output into the ficture transcripts table.

    :param matrix_mtx: str (path) matrix.mtx.gz
    :param barcodes_tsv: str (path) barcodes.tsv.gz
    :param features_tsv: str (path) features.tsv.gz
    :param positions_parquet: str (path) tissue_positions.parquet
    :param scalefactors_json: str (path) scalefactors_json.json
    :param params: dict of config params (feature_type, gene_filter)
    :return: pandas dataframe ["transcript_id", "x", "y", "gene", "count"]
    """
    logging.debug("Reading microns per pixel")
    microns_per_pixel = read_microns_per_pixel(scalefactors_json)
    logging.debug("microns_per_pixel = %s", microns_per_pixel)

    logging.debug("Reading barcodes and features")
    barcodes = read_barcodes(barcodes_tsv)
    features = read_features(features_tsv)

    logging.debug("Reading matrix")
    n_comments, size = count_mtx_header_lines(matrix_mtx)
    logging.debug("Matrix header: %s comment line(s), size %s", n_comments, size)
    mtx = read_matrix(matrix_mtx, n_skip=n_comments + 1)

    # Map 1-based matrix indices to gene symbols and barcodes
    logging.debug("Mapping feature and barcode indices")
    gene_arr = features["gene"].to_numpy()
    mtx["gene"] = gene_arr[mtx["feature_idx"].to_numpy() - 1]
    mtx["barcode"] = barcodes[mtx["barcode_idx"].to_numpy() - 1]

    # Keep only the requested feature type (e.g. "Gene Expression")
    feature_type = params.get("feature_type")
    if feature_type and "feature_type" in features.columns:
        logging.debug("Filtering to feature_type == %s", feature_type)
        ftype_arr = features["feature_type"].to_numpy()
        keep = ftype_arr[mtx["feature_idx"].to_numpy() - 1] == feature_type
        mtx = mtx[keep]

    # Optionally drop junk genes by regex (only if a non-empty pattern is given)
    gene_filter = params.get("gene_filter")
    if gene_filter:
        logging.debug("Dropping genes matching %s", gene_filter)
        mtx = mtx[~mtx["gene"].str.contains(gene_filter)]

    # Attach spatial coordinates and convert pixels to microns
    logging.debug("Joining spatial coordinates")
    positions = read_positions(positions_parquet)
    mtx = mtx.merge(positions, left_on="barcode", right_index=True, how="inner")

    mtx["x"] = mtx["pxl_col_in_fullres"] * microns_per_pixel
    mtx["y"] = mtx["pxl_row_in_fullres"] * microns_per_pixel

    # Synthetic, unique transcript ids (one per non-zero barcode/gene entry)
    mtx = mtx.reset_index(drop=True)
    mtx["transcript_id"] = mtx.index.astype(str)

    return mtx[["transcript_id", "x", "y", "gene", "count"]]


def main(**kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    logging.debug("Running visiumhd_to_tsv.py")

    transcripts_df = visiumhd_to_transcripts(
        matrix_mtx=kwargs["matrix_mtx"],
        barcodes_tsv=kwargs["barcodes_tsv"],
        features_tsv=kwargs["features_tsv"],
        positions_parquet=kwargs["positions_parquet"],
        scalefactors_json=kwargs["scalefactors_json"],
        params=kwargs["params"],
    )

    logging.debug("Writing output transcripts: %s rows", len(transcripts_df))
    transcripts_df.to_pickle(kwargs["out_tsv"])


if __name__ == "__main__":
    main(
        matrix_mtx=snakemake.input.mtx,
        barcodes_tsv=snakemake.input.barcodes,
        features_tsv=snakemake.input.features,
        positions_parquet=snakemake.input.positions,
        scalefactors_json=snakemake.input.scalefactors,
        out_tsv=snakemake.output[0],
        log_file=snakemake.log[0],
        params=snakemake.params.params,
    )
