import logging
import pandas as pd
from sklearn.utils import shuffle

try:
    from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins
    from spers.ficture.workflow.scripts.plot import plot_ficture, plot_factor_score
except ModuleNotFoundError:
    from hex_bin import transcript_to_hex_bins
    from plot import plot_ficture, plot_factor_score


def main(in_scr=None, in_trn=None, out_top=None, out_factors=None,
         log_file=None, params=None, plot=None, platform=None):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    # One per-factor output file per factor; columns are "0".."n-1"
    factor_cols = [str(k) for k in range(len(out_factors))]
    hex_point_scale = params["grid_score"]["hex_width"] ** 2 * 0.9
    # Continuous sequential colour map for the per-factor score plots
    factor_plot = {**plot, "palette": "magma"}

    if platform == "visiumhd":
        # Spot-level scores are already per hex bin with centroid coords + factors
        logging.debug("Reading spot-level hex-bin scores")
        scores_df = pd.read_pickle(in_scr)
        scores_df = scores_df[scores_df["Count"] >= params["plot_filter"]["min_transcripts_per_hex"]]

        logging.debug("Plotting top factor")
        plot_ficture(
            scores_df[["x", "y", "topK"]], out_top,
            point_scale=hex_point_scale, font_scale=0.1, **plot)

        logging.debug("Plotting per-factor scores")
        for factor, out_png in zip(factor_cols, out_factors):
            plot_factor_score(
                scores_df[["x", "y", factor]], out_png, factor=factor,
                point_scale=hex_point_scale, font_scale=0.1, **factor_plot)
        return

    # ---- per-transcript path (xenium / cosmx) ----
    logging.debug("Reading in transcript coords")
    transcripts_df = pd.read_pickle(in_trn)[["transcript_id", "x", "y"]]

    logging.debug("Reading in transcript scores")
    scores_df = pd.read_pickle(in_scr)[["transcript_id", "topK"] + factor_cols]

    logging.debug("Joining transcript classifications")
    transcripts_df = transcripts_df.merge(scores_df, on="transcript_id", how="inner")

    logging.debug("Calculating new hex bins for plotting")
    transcripts_df = transcript_to_hex_bins(transcripts_df, hex_width=params["grid_score"]["hex_width"])

    logging.debug("Filtering low transcript bins")
    hex_counts = transcripts_df.groupby("hex_id").size()
    transcripts_df = transcripts_df[
        transcripts_df["hex_id"].map(hex_counts) >= params["plot_filter"]["min_transcripts_per_hex"]]

    logging.debug("Plotting top factor (down-sampled transcripts)")
    top_df = shuffle(transcripts_df).groupby("hex_id").head(params["plot_filter"]["max_transcripts_per_hex"])
    plot_ficture(top_df, out_top, point_scale=2, font_scale=0.1, **plot)

    logging.debug("Aggregating per-hex factor scores")
    agg_spec = {"x": ("xbin", "first"), "y": ("ybin", "first")}
    agg_spec.update({c: (c, "mean") for c in factor_cols})
    agg = transcripts_df.groupby("hex_id").agg(**agg_spec).reset_index()

    logging.debug("Plotting per-factor scores")
    for factor, out_png in zip(factor_cols, out_factors):
        plot_factor_score(
            agg[["x", "y", factor]], out_png, factor=factor,
            point_scale=hex_point_scale, font_scale=0.1, **factor_plot)


if __name__ == "__main__":
    main(
        in_scr=snakemake.input.scr,
        in_trn=snakemake.input.trn,
        out_top=snakemake.output.top,
        out_factors=list(snakemake.output.factors),
        log_file=snakemake.log[0],
        params=snakemake.params.params,
        plot=snakemake.params.plot,
        platform=snakemake.params.platform
    )
