import seaborn
import matplotlib.pyplot as plt


def plot_ficture(
        plot_df,
        out_png,
        microns_per_inch=100,
        hex_width=16,
        font_scale=0.1,
        point_scale=0.002,
        plot_top_n_factors=12,
        background="black",
        palette="Spectral",
        **params):
    """
    Generate a seaborn plot object for coarse scored transcripts
    :param plot_df: pandas dataframe of fit results ["x", "y", "topK", ...]
    :param out_png: str filepath of output plot
    :param params: config params passed from snakemake rule
    :return: plot object to save
    """

    # Get figure dimensions for output plot
    x_scale = (plot_df["x"].max() - plot_df["x"].min()) / microns_per_inch
    y_scale = (plot_df["y"].max() - plot_df["y"].min()) / microns_per_inch
    plt.figure(figsize=(x_scale, y_scale))

    # Get the font and point scales
    plot_font_scale = max(x_scale, y_scale) * font_scale
    plot_point_scale = (x_scale * y_scale) * hex_width**2 * point_scale

    # Hue order by factor size
    hue_order = list(
        plot_df.groupby("topK").size().sort_values(ascending=False).head(plot_top_n_factors - 1).index)

    # Top n factors
    plot_df.loc[~plot_df["topK"].isin(hue_order), "topK"] = plot_top_n_factors

    # Plot!
    seaborn.set(
        rc={
            "axes.facecolor": background,
            "figure.facecolor": background,
            "axes.grid": False},
        font_scale=plot_font_scale)

    coarse_plot = seaborn.scatterplot(
        plot_df,
        x="x",
        y="y",
        hue="topK",
        s=plot_point_scale,
        palette=palette,
        hue_order=hue_order,
        edgecolor=None)

    coarse_plot.legend([],[], frameon=False)
    coarse_plot.set_xlabel("X (microns)")
    coarse_plot.set_ylabel("Y (microns)")

    plt.savefig(out_png)