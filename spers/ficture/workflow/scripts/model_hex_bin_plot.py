import logging
import pandas as pd

from spers.ficture.workflow.scripts.plot import plot_ficture


def main(plot=None, hex_width=None, **kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)

    logging.debug("Reading in coarse model fit scores")
    fit_df = pd.read_csv(kwargs["in_fit"], sep="\t", compression="gzip", usecols=["hex_id", "topK", "x", "y"])

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
        hex_width=snakemake.params.hex_width
    )