import logging
import pandas as pd

try:
    from spers.ficture.workflow.scripts.plot import plot_ficture
except ModuleNotFoundError:
    from plot import plot_ficture


def main(plot=None, hex_width=None, sample=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)

    logging.debug("Reading in coarse model fit scores")
    fit_df = pd.read_csv(kwargs["in_fit"], sep="\t", compression="gzip", usecols=["sample", "hex_id", "topK", "x", "y"])

    # The fit is a combined (joint) table; plot just this sample's hex bins
    logging.debug("Selecting sample %s", sample)
    fit_df = fit_df[fit_df["sample"] == sample]

    logging.debug("Plotting model hex IDs")
    plot_ficture(
        fit_df,
        kwargs["out_png"],
        font_scale=0.1,
        point_scale=hex_width**2 * 0.9,
        **plot
    )


if __name__ == "__main__":
    main(
        in_fit=snakemake.input.fit,
        out_png=snakemake.output.png,
        log_file=snakemake.log[0],
        threads=snakemake.threads,
        plot=snakemake.params.plot,
        hex_width=snakemake.params.hex_width,
        sample=snakemake.params.sample
    )