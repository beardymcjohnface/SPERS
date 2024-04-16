import logging, gzip
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


def batch_dimensions(df, **params):
    """
    Calculate the x and y batch spacing to approximately accommodate the batch size

    :param df: Pandas DF of transcripts ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    :param params: dict of config params ["batch_size", "batch_buff"]
    :return: params dict updated with "x_step", "y_step", "x_max", "y_max" ...
    """

    # Finding x and y min and max values
    params["x_min"], params["y_min"] = df[["xbin", "ybin"]].values.min(axis=0)
    params["x_max"], params["y_max"] = df[["xbin", "ybin"]].values.max(axis=0)

    # Find x and y ranges
    params["xrange"] = params["x_max"] - params["x_min"]
    params["yrange"] = params["y_max"] - params["y_min"]

    # Adjust x and y batch sizes to fit data ranges
    params["x_batch"] = round(params["xrange"] / (params["batch_size"] - params["batch_buff"]))
    params["y_batch"] = round(params["yrange"] / (params["batch_size"] - params["batch_buff"]))

    # Calculate x and y batch steps
    params["x_step"] = ( params["xrange"] / params["x_batch"] ) - params["batch_buff"]
    params["y_step"] = ( params["yrange"] / params["y_batch"] ) - params["batch_buff"]

    return params


def minibatch_transcripts(df, file, **params):
    """
    Create batches for transcripts based on batch dimensions

    :param df: transcripts pandas dataframe ["hex_id", "transcript_id", "xbin", "ybin", "gene"]
    :param file: output file path (str)
    :param params: dict config params
    :return: None
    """
    
    # Write header to output file
    header = ["minibatch_id", "hex_id"]
    with gzip.open(file, "wt") as wf:
        wf.write("\t".join(header) + "\n")

    # Generate random integers for batch IDs
    rand_ints = set(np.random.randint(1000000000, size=1000000))

    # Iterate x and y steps and batch on the fly
    x_start = params["x_min"]
    while x_start + params["x_step"] + params["batch_buff"] <= params["x_max"]:
        y_start = params["y_min"]
        while y_start + params["y_step"] + params["batch_buff"] <= params["y_max"]:

            # Get dataframe slice for transcripts within batch boundaries
            df_slice = df[(df.xbin >= x_start) & (df.xbin <= x_start + params["x_step"] + params["batch_buff"]) &
                          (df.ybin >= y_start) & (df.ybin <= y_start + params["y_step"] + params["batch_buff"])]

            # Add random index for the batch
            df_slice["minibatch_id"] = rand_ints.pop()

            # Append the slice to the output file
            df_slice[header].drop_duplicates().to_csv(file, mode="a", sep="\t", compression="gzip", header=False, float_format="%.2f")

            y_start += params["y_step"]
        x_start += params["x_step"]



def main(infile=None, outfile=None, log_file=None, params=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading input transcripts: " + str(infile))
    transcripts_df = read_input(infile)

    logging.debug("Calculating minibatch dimensions")
    params = batch_dimensions(transcripts_df, **params)

    logging.debug("Batching and writing transcripts")
    minibatch_transcripts(transcripts_df, outfile, **params)


if __name__ == "__main__":
    main(
        infile=snakemake.input[0],
        outfile=snakemake.output[0],
        log_file=snakemake.log[0],
        params=snakemake.params.params
    )