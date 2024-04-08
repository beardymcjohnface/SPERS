import logging
import pandas as pd


def read_transcripts(in_csv, params):
    """
    Read in the CosMX CSV-formatted transcripts file and convert to TSV

    :param in_csv: str (path) CSV file optionally compressed with gzip
    :param params: params for parsing and processing
    :return: pandas dataframe
    """
    # Read in file
    if in_csv.endswith(".gz"):
        df = pd.read_csv(in_csv, compression="gzip")
    else:
        df = pd.read_csv(in_csv)

    # Remove unnecessary columns and rename
    df = df[[params["headers"]["X"], params["headers"]["Y"], params["headers"]["gene"]]]
    df.columns = ["X", "Y", "gene"]

    # Filter junk rows
    df = df[~df["gene"].str.contains("Control")]
    df = df[~df["gene"].str.contains("Neg")]
    df = df[~df["gene"].str.contains("Unassigned")]

    # Group and sort - add counts
    df = df.groupby(["X", "Y", "gene"]).size().reset_index(name="Count")
    df = df.sort_values(by="Y")

    # Rescale XY offset to zero
    df["Y"] -= df["Y"].min()
    df["X"] -= df["X"].min()

    # Rescale into um
    df["X"] /= params["mu_scale"]
    df["Y"] /= params["mu_scale"]

    return df


def coordinate_minmax(df):
    """
    Collect the Min/Max X/Y values of the converted transcripts file

    :param df: pandas dataframe from read_transcripts func
    :return: pandas dataframe of x/y min/max for writing
    """
    xmin = df["X"].min()
    xmax = df["X"].max()
    ymin = df["Y"].min()
    ymax = df["Y"].max()

    minmax_df = {"xmin": [xmin], "xmax": [xmax], "ymin": [ymin], "ymax": [ymax]}
    minmax_df = pd.DataFrame(minmax_df).transpose()

    return minmax_df


def count_features(df):
    """
    Sum the counts for all genes

    :param df: pandas dataframe from read_transcripts func
    :return:
    """
    df_feature = df.drop(["X", "Y"], axis=1).copy()
    df_feature = df_feature.groupby("gene")["Count"].sum().reset_index()

    return df_feature


def main(**kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    logging.debug("Running transcript_csv_to_tsv.py")

    logging.debug("Reading and converting CSV")
    transcripts_df = read_transcripts(kwargs["input_csv"], kwargs["params"])

    logging.debug("Collecting min/max values")
    minmax_coords = coordinate_minmax(transcripts_df)

    logging.debug("Collecting gene summaries")
    feature_counts = count_features(transcripts_df)

    logging.debug("Writing output transcripts")
    transcripts_df.to_csv(kwargs["out_tsv"], sep="\t", compression="gzip", index=False, float_format="%.2f")

    logging.debug("Writing output gene counts")
    feature_counts.to_csv(kwargs["out_feat"], sep="\t", compression="gzip", index=False)

    logging.debug("Writing output min/max coords")
    minmax_coords.to_csv(kwargs["out_minmax"], sep="\t", index=True, header=False, float_format="%.2f")


if __name__ == "__main__":
    main(
        input_csv=snakemake.input[0],
        out_tsv=snakemake.output.tsv,
        out_feat=snakemake.output.feat,
        out_minmax=snakemake.output.minmax,
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )
