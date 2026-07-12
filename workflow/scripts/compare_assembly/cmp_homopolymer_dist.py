import sys
import argparse

import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from scipy.stats import normaltest, ks_2samp, false_discovery_control

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["text.usetex"] = False

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


def compare_homopolymer_bins(df_nt_bins: pl.DataFrame) -> dict[tuple[str, str], float]:
    dfs_nt_label_bins = df_nt_bins.partition_by(["label"], as_dict=True)
    for lbl, df_vals in dfs_nt_label_bins.items():
        lbl = lbl[0]
        nres = normaltest(df_vals["len"])
        assert nres.pvalue < 0.05, "One of labels is normally distributed"
        # dim = df_vals.filter(pl.col("len") == pl.col("len").mean().cast(pl.Int64)).shape
        # print(f"{lbl} mean length: {df_vals['len'].mean()} ({dim[0]})", file=sys.stderr)
        # print(f"{lbl}: {nres}", file=sys.stderr)

    # Use Komologorov-smirnov test and post-hoc analysis to determine if one assembly has more homopolymer errors than another assembly.
    # We use this instead of our original kruskal wallis as the distribution is compared rather than just the median.
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
    pvals = {}
    for lbl, df in dfs_nt_label_bins.items():
        lbl = lbl[0]
        for lbl2, df2 in dfs_nt_label_bins.items():
            lbl2 = lbl2[0]
            res = ks_2samp(df["len"], df2["len"], alternative="greater")
            pvals[(lbl, lbl2)] = res.pvalue

    # Requires list for some reason
    adj_pvals = false_discovery_control(list(pvals.values()))
    # print(pvals, file=sys.stderr)
    # print(adj_pvals, file=sys.stderr)
    return {
        lbl: adj_pval
        for (lbl, pval), adj_pval in zip(pvals.items(), adj_pvals, strict=True)
    }


# https://rowannicholls.github.io/python/graphs/ax_based/boxplots_significance.html
def draw_signif_brackets(
    ax: Axes, x1: int, x2: int, level: float, ylim: tuple[float, float], p: float
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
    ax.text(
        (x1 + x2) * 0.5,
        text_height,
        sig_symbol,
        ha="center",
        va="bottom",
        c="k",
        fontsize=14,
    )


def draw_plot_all(
    df_bins: pl.DataFrame,
    labels: list[str],
    label_idxs: dict[str, int],
    label_colors: dict[str, str],
    figsize: tuple[int, int] = (8, 8),
):
    ylims = [(350, 400), (0, 100)]
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=figsize,
        layout="constrained",
        height_ratios=[(ylim[1] - ylim[0]) / 100.0 for ylim in ylims],
    )
    ax: Axes
    ax1: Axes
    ax2: Axes
    alpha = 0.05
    ax1 = axes[0]
    ax2 = axes[1]
    res_map = compare_homopolymer_bins(df_bins)

    for ax, ylim in (
        (ax1, ylims[0]),
        (ax2, ylims[1]),
    ):
        sns.violinplot(
            data=df_bins,
            x="label",
            y="len",
            hue="label",
            density_norm="count",
            inner="quart",
            legend=None,
            palette=label_colors.values(),
            hue_order=labels,
            order=labels,
            alpha=0.5,
            ax=ax,
            cut=0,
        )
        sns.stripplot(
            data=df_bins,
            x="label",
            y="len",
            hue="label",
            jitter=0.001,
            palette=label_colors.values(),
            linewidth=1,
            hue_order=labels,
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
        level += 2

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
    ax1.tick_params(axis="both", which="major", labelsize=14)
    ax2.tick_params(axis="both", which="major", labelsize=14)
    ax1.set_xticks([], [])
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    add_slashes(ax1, ax2)

    for lbl in ax.xaxis.get_majorticklabels():
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")
        lbl.set_fontsize(14)
        lbl.set_color(label_colors.get(lbl.get_text(), "white"))

    fig.supylabel("Error homopolymer length (bp)", fontsize=14)
    return fig


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
    ap.add_argument("-c", "--colors", type=str, nargs="+")
    ap.add_argument(
        "-o",
        "--output_prefix",
        default="all",
        help="Plot prefix for all overlapped.",
    )

    args = ap.parse_args()

    labels = args.labels
    label_colors = dict(zip(labels, args.colors, strict=True))
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
        res_map = compare_homopolymer_bins(df_nt_bins)

        ax1.set_title(nt, fontsize=18)

        for ax, ylim in (
            (ax1, ylims[0]),
            (ax2, ylims[1]),
        ):
            sns.violinplot(
                data=df_nt_bins,
                x="label",
                y="len",
                color=color,
                inner="quart",
                legend=None,
                density_norm="count",
                order=labels,
                hue_order=labels,
                alpha=0.5,
                cut=0,
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
                hue_order=labels,
                legend=None,
                alpha=0.5,
                ax=ax,
            )
            ax.set_ylim(ylim)
            ax.set_ylabel(None)

        level = 0
        ymin, ymax = ax1.get_ylim()
        # Draw brackets within this yaxis limit
        ylim = (ymin, ymax - 12)
        for (label_1, label_2), p in res_map.items():
            if p > alpha:
                print(
                    f"No significant difference between {nt} {label_1} and {label_2} with {p=}",
                    file=sys.stderr,
                )
                continue
            x1 = label_idxs[label_1]
            x2 = label_idxs[label_2]
            # The alternative is that F(x) > G(x) for at least one x.
            print(
                f"(*) Significant difference for {nt} where {label_2} has more error homopolymers at some homopolymer length than {label_1} with {p=}",
                file=sys.stderr,
            )
            draw_signif_brackets(ax=ax1, x1=x1, x2=x2, level=level, ylim=ylim, p=p)
            level += 2

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
        ax1.tick_params(axis="both", which="major", labelsize=14)
        ax2.tick_params(axis="both", which="major", labelsize=14)
        ax1.set_xticks([], [])
        ax1.set_xlabel(None)
        ax2.set_xlabel(None)
        add_slashes(ax1, ax2)

        for lbl in ax.xaxis.get_majorticklabels():
            lbl.set_rotation(45)
            lbl.set_fontsize(14)
            lbl.set_horizontalalignment("right")
            lbl.set_rotation_mode("anchor")
            lbl.set_color(label_colors.get(lbl.get_text(), "white"))

    # fig.legend(
    #     handles=[
    #         Patch(facecolor=color, edgecolor="black", label=lbl, alpha=0.5)
    #         for lbl, color in NT_COLORS.items()
    #     ],
    #     loc="center left",
    #     bbox_to_anchor=(1, 0.9),
    #     fontsize=14,
    #     **LEGEND_KWARGS,
    # )
    fig.supylabel("Error homopolymer length (bp)", fontsize=14)
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=300)

    fig_all = draw_plot_all(df_bins, labels, label_idxs, label_colors)
    fig_all.savefig(f"{args.output_prefix}_all.png", bbox_inches="tight", dpi=300)
    fig_all.savefig(f"{args.output_prefix}_all.pdf", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
