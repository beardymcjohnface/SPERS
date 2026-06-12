import logging
import pandas as pd

try:
    from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
except ModuleNotFoundError:
    from hex_bin import transcript_to_hex_bins


def main(in_tsv=None, out_pkl=None, hex_width=None, sample=None, log_file=None):
    """
    Hex-bin a single sample and dump a compact, pre-aggregated representation so the
    joint training step never has to hold every sample's raw transcripts at once.

    Writes a pickle of a dict:
      - "counts": long dataframe ["hex_id", "gene", "count"] (hex x gene counts)
      - "meta":   one row per bin ["hex_id", "sample", "xbin", "ybin"]
    Hex ids are namespaced by sample ("<sample>:<local>") so they never collide
    when pooled across samples.
    """
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading transcripts for sample %s", sample)
    df = pd.read_pickle(in_tsv)

    logging.debug("Hex binning")
    df = transcript_to_hex_bins(df, x_offset=0, y_offset=0, hex_width=hex_width)
    df["hex_id"] = str(sample) + ":" + df["hex_id"].astype(str)

    logging.debug("Aggregating hex x gene counts")
    if "count" in df.columns:
        counts = df.groupby(["hex_id", "gene"])["count"].sum().reset_index()
    else:
        counts = df.groupby(["hex_id", "gene"]).size().reset_index()
    counts.columns = ["hex_id", "gene", "count"]

    meta = df[["hex_id", "xbin", "ybin"]].drop_duplicates("hex_id")
    meta["sample"] = str(sample)

    logging.debug("Writing %d bins / %d (hex,gene) entries", len(meta), len(counts))
    pd.to_pickle({"counts": counts, "meta": meta}, out_pkl)


if __name__ == "__main__":
    main(
        in_tsv=snakemake.input[0],
        out_pkl=snakemake.output[0],
        hex_width=snakemake.params.hex_width,
        sample=snakemake.params.sample,
        log_file=snakemake.log[0],
    )
