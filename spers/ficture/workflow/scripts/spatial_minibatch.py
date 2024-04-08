import logging, gzip
import numpy as np
import pandas as pd
import sys


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
    Calculate the X and Y batch spacing to approximately accommodate the batch size

    :param df: Pandas DF of transcripts ["X", "Y", "gene", "Count"]
    :param params: dict of config params ["batch_size", "batch_buff"]
    :return: params dict updated with "x_step", "y_step", "x_max", "y_max" ...
    """

    # Finding X and Y min and max values
    params["x_min"], params["y_min"] = df[["X", "Y"]].values.min(axis=0)
    params["x_max"], params["y_max"] = df[["X", "Y"]].values.max(axis=0)

    # Find X and Y ranges
    params["xrange"] = params["x_max"] - params["x_min"]
    params["yrange"] = params["y_max"] - params["y_min"]

    # Adjust X and Y batch sizes to fit data ranges
    params["x_batch"] = round(params["xrange"] / (params["batch_size"] - params["batch_buff"]))
    params["y_batch"] = round(params["yrange"] / (params["batch_size"] - params["batch_buff"]))

    # Calculate X and Y batch steps
    params["x_step"] = ( params["xrange"] / params["x_batch"] ) - params["batch_buff"]
    params["y_step"] = ( params["yrange"] / params["y_batch"] ) - params["batch_buff"]

    return params


def minibatch_transcripts(df, file, **params):
    """
    Create batches for transcripts based on batch dimensions

    :param df: transcripts pandas dataframe ["X", "Y", "gene", "Count"]
    :param file: output file path (str)
    :param params: dict config params
    :return: None
    """
    
    # Write header to output file
    header = ["random_index", "X", "Y", "Count"]
    with gzip.open(file, "wt") as wf:
        wf.write("\t".join(header) + "\n")

    # Generate random integers for batch IDs
    rand_ints = set(np.random.randint(1000000000, size=1000000))

    # Iterate X and Y steps and batch on the fly
    x_start = params["x_min"]
    while x_start + params["x_step"] + params["batch_buff"] <= params["x_max"]:
        y_start = params["y_min"]
        while y_start + params["y_step"] + params["batch_buff"] <= params["y_max"]:

            # Get dataframe slice for transcripts within batch boundaries
            df_slice = df[(df.X >= x_start) & (df.X <= x_start + params["x_step"] + params["batch_buff"]) &
                          (df.Y >= y_start) & (df.Y <= y_start + params["y_step"] + params["batch_buff"])]

            # Add random index for the batch
            df_slice["random_index"] = rand_ints.pop()

            # Append the slice to the output file
            df_slice[header].to_csv(file, mode="a", sep="\t", compression="gzip", header=False, float_format="%.2f")

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