import sys
import argparse

import polars as pl
import seaborn as sns
import scikit_posthocs as sp
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from scipy.stats import normaltest

NT_COLORS = {
    "A": "#009600",
    "C": "#0000FF",
    "G": "#D17105",
    "T": "#FF0000",
}
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)


# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
def add_slashes(ax1: Axes, ax2: Axes):
    d = 0.5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(
        marker=[(-1, -d), (1, d)],
        markersize=12,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )
    ax1.plot(0, 0, transform=ax1.transAxes, **kwargs)
    ax1.get_yticklines()[0].set_visible(False)
    ax2.plot(0, 1, transform=ax2.transAxes, **kwargs)

    ymax = ax2.get_yticklines()[-1].get_ydata()[0]
    for ytickline in ax2.get_yticklines():
        y = ytickline.get_ydata()[0]
        if y == ymax:
            ytickline.set_visible(False)


def dunn_test(df_nt_bins: pl.DataFrame) -> dict[tuple[str, str], float]:
    dfs_nt_label_bins = df_nt_bins.partition_by(["label"], as_dict=True)
    for lbl, df_vals in dfs_nt_label_bins.items():
        lbl = lbl[0]
        nres = normaltest(df_vals["len"])
        assert nres.pvalue < 0.05, "One of labels is normally distributed"
        print(f"{lbl} mean length: {df_vals['len'].mean()}", file=sys.stderr)
        print(f"{lbl}: {nres}", file=sys.stderr)

    # Use Kruskal-Wallis test and post-hoc Dunn to determine if any difference between medians of labels.
    # https://stats.stackexchange.com/a/95270
    # https://www.statology.org/dunns-test-python/
    # H_0 - There is no significant difference in the median between a pair of labels.
    res = sp.posthoc_dunn(
        [df["len"] for df in dfs_nt_label_bins.values()], p_adjust="bonferroni"
    )
    print(res, file=sys.stderr)
    labels = tuple(dfs_nt_label_bins.keys())
    res.index = labels
    res.columns = labels
    return {
        tuple(sorted([label_1[0], label_2[0]])): res.loc[label_1, label_2]
        for label_1 in res.index
        for label_2 in res.columns
        if label_1 != label_2
    }


# https://rowannicholls.github.io/python/graphs/ax_based/boxplots_significance.html
def draw_signif_brackets(
    ax: Axes, x1: int, x2: int, level: int, ylim: tuple[float, float], p: float
):
    # What level is this bar among the bars above the plot?
    # Plot the bar
    y_range = ylim[1] - ylim[0]
    bar_height = (y_range * 0.09 * level) + ylim[1]
    bar_tips = bar_height - (y_range * 0.02)
    ax.plot(
        [x1, x1, x2, x2],
        [bar_tips, bar_height, bar_height, bar_tips],
        lw=1,
        c="k",
        clip_on=False,
    )
    # Significance level
    if p < 0.001:
        sig_symbol = "***"
    elif p < 0.01:
        sig_symbol = "**"
    elif p < 0.05:
        sig_symbol = "*"
    text_height = bar_height + (y_range * 0.01)
    ax.text((x1 + x2) * 0.5, text_height, sig_symbol, ha="center", va="bottom", c="k")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--bins",
        nargs="+",
        type=argparse.FileType("rb"),
        help="Homopolymer bins TSV file.",
    )
    ap.add_argument(
        "--labels",
        nargs="+",
        type=str,
        help="Homopolymer BED4 file labels",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="all.png",
        help="Plot for all overlapped.",
    )

    args = ap.parse_args()

    labels = args.labels
    label_idxs = {label: i for i, label in enumerate(labels)}
    dfs_bins: list[pl.DataFrame] = []
    for label, file in zip(labels, args.bins, strict=True):
        df_bins = pl.read_csv(
            file,
            has_header=False,
            separator="\t",
            new_columns=["nt", "len", "count"],
        ).with_columns(label=pl.lit(label))

        df_bins = pl.DataFrame(
            [
                {"nt": nt, "label": label, "len": ln}
                for nt, ln, count, label in df_bins.iter_rows()
                for _ in range(count)
                # Ignore non-overlapping longdust or NucFlag regions.
                if ln != 0
            ],
            orient="row",
        )

        dfs_bins.append(df_bins)

    df_bins = pl.concat(dfs_bins)

    ylims = [(350, 400), (0, 100)]

    fig, axes = plt.subplots(
        nrows=2,
        ncols=len(NT_COLORS.keys()),
        figsize=(16, 8),
        layout="constrained",
        height_ratios=[(ylim[1] - ylim[0]) / 100.0 for ylim in ylims],
    )
    ax: Axes
    ax1: Axes
    ax2: Axes
    alpha = 0.05
    for i, (nt, color) in enumerate(NT_COLORS.items()):
        ax1 = axes[0, i]
        ax2 = axes[1, i]
        df_nt_bins = df_bins.filter(pl.col("nt") == nt)
        res_map = dunn_test(df_nt_bins)

        for ax, ylim in (
            (ax1, ylims[0]),
            (ax2, ylims[1]),
        ):
            sns.violinplot(
                data=df_nt_bins,
                x="label",
                y="len",
                color=color,
                split=True,
                gap=0.5,
                inner="quart",
                legend=None,
                order=labels,
                alpha=0.5,
                ax=ax,
            )
            sns.stripplot(
                data=df_nt_bins,
                x="label",
                y="len",
                color=color,
                dodge=True,
                jitter=0.001,
                linewidth=1,
                order=labels,
                legend=None,
                alpha=0.5,
                ax=ax,
            )
            ax.set_ylim(ylim)
            ax.set_ylabel(None)

        level = 0
        for (label_1, label_2), p in res_map.items():
            if p > alpha:
                continue
            x1 = label_idxs[label_1]
            x2 = label_idxs[label_2]
            draw_signif_brackets(
                ax=ax1, x1=x1, x2=x2, level=level, ylim=ax1.get_ylim(), p=p
            )
            level += 1

        for spine in ("top", "right", "bottom"):
            ax1.spines[spine].set_visible(False)

        for spine in ("top", "right"):
            ax2.spines[spine].set_visible(False)

        ax1.tick_params(
            axis="x",
            which="both",
            bottom=False,
        )
        ax2.tick_params(axis="x", labelrotation=45)
        ax1.set_xticks([], [])
        ax1.set_xlabel(None)
        ax2.set_xlabel(None)
        add_slashes(ax1, ax2)

    fig.legend(
        handles=[
            Patch(color=color, label=lbl, alpha=0.5) for lbl, color in NT_COLORS.items()
        ],
        loc="center left",
        bbox_to_anchor=(1, 0.9),
        **LEGEND_KWARGS,
    )
    fig.supylabel("Homopolymer length (bp)")
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
