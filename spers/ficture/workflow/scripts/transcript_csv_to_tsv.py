import logging
import pandas as pd


def read_transcripts(in_csv, params):
    """
    Read in the CosMX CSV-formatted transcripts file and convert to TSV

    :param in_csv: str (path) CSV file optionally compressed with gzip
    :param params: params for parsing and processing
    :return: pandas dataframe ["transcript_id", "x", "y", "gene"]
    """
    # Read in file
    if in_csv.endswith(".gz"):
        df = pd.read_csv(in_csv, compression="gzip")
    else:
        df = pd.read_csv(in_csv)

    # Remove unnecessary columns and rename
    df = df[list(params["headers"].values())]
    df.columns = list(params["headers"].keys())

    # Add transcript index if missing
    if "transcript_id" not in params["headers"]:
        df["transcript_id"] = df.index
        df["transcript_id"] = df["transcript_id"].astype(str).str.zfill(15)

    # Filter junk rows
    df = df[~df["gene"].str.contains(params["gene_filter"])]

    # # Group and sort - add counts
    # df = df.groupby(["x", "y", "gene"]).size().reset_index(name="Count")
    # df = df.sort_values(by="y")

    # Rescale Xy offset to zero
    # df["y"] -= df["y"].min()
    # df["x"] -= df["x"].min()

    # Rescale into um
    df["x"] /= params["mu_scale"]
    df["y"] /= params["mu_scale"]

    return df


def main(**kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    logging.debug("Running transcript_csv_to_tsv.py")

    logging.debug("Reading and converting CSV")
    transcripts_df = read_transcripts(kwargs["input_csv"], kwargs["params"])

    logging.debug("Writing output transcripts")
    transcripts_df.to_csv(kwargs["out_tsv"], sep="\t", compression="gzip", index=False, float_format="%.2f")



if __name__ == "__main__":
    main(
        input_csv=snakemake.input[0],
        out_tsv=snakemake.output[0],
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )
