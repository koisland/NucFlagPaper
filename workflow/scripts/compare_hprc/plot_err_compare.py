import argparse
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from typing import Literal
from matplotlib.patches import Patch
from matplotlib.colors import rgb2hex
from matplotlib.axes import Axes

COLS_CALLS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "tstart",
    "tend",
    "item_rgb",
)
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    title=None,
    fancybox=False,
    frameon=False,
)
DIPLOID_GENOME_SIZE = 6_200_000_000
AFR_POP_LABELS = {"ACB", "GWD", "ESN", "MKK", "LWK", "ASL", "YRI", "MSL", "ASW"}
POP_GROUP_COLORS = {
    "AFR": "orange",
    "Non-AFR": "gray",
}
LARGE_ERRORS = (
    "false_dup",
    "misjoin",
    "collapse",
    "scaffold",
    "softclip",
)
SMALL_ERRORS = (
    "insertion",
    "deletion",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "other_repeat",
    "mismatch",
    "low_quality",
    "het_or_mismap",
)
plt.rcParams["font.family"] = "Arial"


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


def get_ylim_w_buffer(
    df: pl.DataFrame, col: str, buffer: float = 0.1
) -> tuple[float, float]:
    ylim_min, ylim_max = df[col].min(), df[col].max()
    # Add 10% of max as buffer
    ylim_buffer = buffer * max(abs(ylim_min), abs(ylim_max))
    ylim_min += -ylim_buffer
    ylim_max += ylim_buffer
    return ylim_min, ylim_max


def draw_combined_err_diff(
    df: pl.DataFrame,
    output_prefix: str,
    name_colors: dict[str, str],
    include_all: bool = False,
    by: Literal["count", "length"] = "count",
):
    if by == "count":
        col = "diff_count"
        ylabel = "Difference in number of calls between R1 and R2"
    else:
        col = "diff"
        ylabel = "Difference in length of calls between R1 and R2 (Mbp)"

    df_filtered = df
    if include_all:
        x_order = ["all", *LARGE_ERRORS, *SMALL_ERRORS]
        name_colors = {"all": "#808080", **name_colors}
    else:
        x_order = [*LARGE_ERRORS, *SMALL_ERRORS]
        df_filtered = df.filter(pl.col("name").ne(pl.lit("all")))

    # Split axes if all since will be overwhelmingly large
    split_axes = include_all and by == "count"

    if split_axes:
        name_fig, axes_fig = plt.subplots(
            layout="constrained", figsize=(12, 10), nrows=2, height_ratios=[0.5, 0.33]
        )
        name_fig.supylabel(ylabel, fontsize=18)
        # Add slashes and separate ylimits
        add_slashes(axes_fig[0], axes_fig[1])
        ylim_non_all = get_ylim_w_buffer(
            df_filtered.filter(pl.col("name").ne(pl.lit("all"))), col, buffer=0.25
        )
        min_val_all = (
            df_filtered.filter(pl.col("name").eq(pl.lit("all")))
            .max()
            .get_column(col)
            .first()
        )
        ylim_all_length = ylim_non_all[1] - ylim_non_all[0]
        ylim_all = (
            min_val_all - ylim_all_length * 0.25,
            min_val_all + ylim_all_length * 0.5,
        )
        set_ylim = {0: ylim_non_all, 1: ylim_all}
    else:
        name_fig, ax_fig = plt.subplots(layout="constrained", figsize=(12, 10))
        axes_fig = [ax_fig]
        set_ylim = None

    for i, ax_fig in enumerate(axes_fig):
        ax_fig: Axes
        sns.barplot(
            data=df_filtered,
            x="name",
            y=col,
            hue="name",
            palette=name_colors,
            order=x_order,
            legend=None,
            ax=ax_fig,
        )
        for cont in ax_fig.containers:
            kwargs = {}
            if by == "length":
                kwargs = {"fmt": lambda length: f"{length / 1_000_000:.1f}"}
            else:
                kwargs = {"fmt": lambda count: f"{count:,.0f}"}

            ax_fig.bar_label(
                cont,
                fontsize=12,
                label_type="edge",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
                padding=3,
                **kwargs,
            )

        if set_ylim and set_ylim.get(i):
            ylim = set_ylim[i]
        else:
            ylim = get_ylim_w_buffer(df_filtered, col)

        ax_fig.set_ylim(ylim)
        ax_fig.tick_params(axis="both", which="major", labelsize=18)
        ax_fig.margins(x=0.01)
        ax_fig.set_xlabel(None)
        # Don't duplicate label if split axes
        if not split_axes:
            ax_fig.set_ylabel(ylabel, fontsize=18)
        else:
            ax_fig.set_ylabel(None)

        # Only label bottom
        if i == len(axes_fig) - 1:
            for label in ax_fig.get_xticklabels():
                label.set_color(name_colors.get(label.get_text(), "gray"))
                label.set_path_effects(
                    [pe.withStroke(linewidth=0.2, foreground="black")]
                )
                label.set_rotation(45)
                label.set_horizontalalignment("right")
                label.set_rotation_mode("anchor")
            spines = ["top", "right"]
        else:
            spines = ["top", "bottom", "right"]
            ax_fig.tick_params(
                axis="both",
                bottom=False,
                labelbottom=False,
            )

        for spine in spines:
            ax_fig.spines[spine].set_visible(False)

        if by == "length":
            ax_fig.yaxis.set_major_formatter(lambda x, _: f"{x / 1_000_000:.1f}")
        else:
            ax_fig.yaxis.set_major_formatter(lambda x, _: f"{x:,.0f}")

    name_fig.savefig(f"{output_prefix}_combined.png", bbox_inches="tight", dpi=300)
    name_fig.savefig(f"{output_prefix}_combined.pdf", bbox_inches="tight", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-a",
        "--input_a",
        type=str,
        help="Input calls for group a. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument(
        "-b",
        "--input_b",
        type=str,
        help="Input calls for group b. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument("-m", "--metadata", type=str, help="Input sample metadata TSV.")
    ap.add_argument("-o", "--output_prefix", type=str, default="out")

    args = ap.parse_args()

    df_a = (
        pl.read_csv(
            args.input_a,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_CALLS,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("a"), length=pl.col("end") - pl.col("start"))
        .filter(pl.col("name").ne(pl.lit("correct")))
    )
    df_b = (
        pl.read_csv(
            args.input_b,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_CALLS,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("b"), length=pl.col("end") - pl.col("start"))
        .filter(pl.col("name").ne(pl.lit("correct")))
    )
    df_metadata = pl.read_csv(
        args.metadata, separator="\t", has_header=True
    ).with_columns(
        pop_group=pl.when(pl.col("Population Abbreviation").is_in(AFR_POP_LABELS))
        .then(pl.lit("AFR"))
        .otherwise(pl.lit("Non-AFR"))
    )
    # afr_samples: set[Any] = set(
    #     df_metadata.filter(pl.col("pop_group").eq(pl.lit("AFR")))["Sample ID"]
    # )

    df_all = pl.concat((df_a, df_b))
    df_grp_no_hap = df_all.group_by(["sample", "label", "name"]).agg(
        length=pl.col("length").sum(), count=pl.col("length").count().cast(pl.Int32)
    )
    df_grp_no_hap_wide = (
        df_grp_no_hap.pivot(on="name", index=("sample", "label"), values="length")
        .sort(by=("sample", "label"))
        .fill_null(0)
    )

    unique_samples = df_grp_no_hap_wide["sample"].unique().sort()
    name_colors = {
        name: rgb2hex(tuple(int(c) / 255 for c in color.split(",")))
        for name, color in df_all.select("name", "item_rgb").iter_rows()
    }
    label_samples = np.arange(len(unique_samples))
    bottoms = {
        "a": np.zeros(len(unique_samples)),
        "b": np.zeros(len(unique_samples)),
    }

    fig, ax = plt.subplots(layout="constrained", figsize=(24, 6))
    ax: Axes
    width = 0.3
    multiplier_adj = 1.1

    for call_name, call_color in name_colors.items():
        multiplier = 0
        for label in bottoms.keys():
            bottom = bottoms[label]
            offset = width * multiplier
            lengths = df_grp_no_hap_wide.filter(pl.col("label").eq(label))[
                call_name
            ].to_numpy()
            rects = ax.bar(
                label_samples + offset,
                lengths,
                width,
                label=call_name,
                bottom=bottom,
                facecolor=call_color,
                edgecolor="black",
            )
            # Don't plot small regions which you wouldn't see anyways.
            if lengths.mean() > 10_000:
                ax.bar_label(
                    rects,
                    fmt=lambda length: f"{length / 1_000_000:.1f}",
                    fontsize=8,
                    label_type="edge",
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
                    padding=3,
                )

            # Label colors should be red/blue for r1/r2 based on palette
            multiplier += multiplier_adj
            bottom += lengths

    # Minimize spacing between xtick and spine
    ax.margins(x=0.01)

    # Backwards since last one is first one in
    ax.legend(
        labels=reversed(name_colors.keys()),
        handles=[
            Patch(edgecolor="black", facecolor=color)
            for color in reversed(name_colors.values())
        ],
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        **LEGEND_KWARGS,
    )
    # Nudge xtick slightly so centered
    # Set sample names
    ax.set_xticks(label_samples + 0.15, unique_samples)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Color by population group and rotate
    for lbl in ax.xaxis.get_majorticklabels():
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")
        # color = POP_GROUP_COLORS["AFR" if lbl.get_text() in afr_samples else "Non-AFR"]
        # lbl.set_color(color)

    # Just eyeballing
    ax.yaxis.set_major_formatter(lambda x, _: f"{x / 1_000_000:.1f}")
    ax.set_ylabel("Cumulative length (Mbp)")
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=300)

    # Agg difference in called bases by sample and sample/name
    df_cmp_ab = (
        df_grp_no_hap.group_by(["sample", "name", "label"])
        .agg(
            length=pl.col("length").sum(),
            count=pl.col("count").sum(),
        )
        .pivot(on="label", index=["sample", "name"])
        .select("sample", "name", "length_a", "length_b", "count_a", "count_b")
        .with_columns(
            diff=pl.col("length_b") - pl.col("length_a"),
            diff_count=pl.col("count_b") - pl.col("count_a"),
        )
        .join(df_metadata, left_on="sample", right_on="Sample ID", how="left")
    )
    df_cmp_name_ab = (
        df_cmp_ab.group_by(["name"])
        .agg(
            a=pl.col("length_a").sum(),
            a_count=pl.col("count_a").sum(),
            b=pl.col("length_b").sum(),
            b_count=pl.col("count_b").sum(),
        )
        .with_columns(
            diff=pl.col("b") - pl.col("a"),
            diff_count=pl.col("b_count") - pl.col("a_count"),
        )
        .sort(by="name")
    )
    df_cmp_name_ab = pl.concat(
        [
            df_cmp_name_ab,
            df_cmp_name_ab.group_by(None)
            .agg(
                pl.col("a").sum(),
                pl.col("a_count").sum(),
                pl.col("b").sum(),
                pl.col("b_count").sum(),
                pl.col("diff").sum(),
                pl.col("diff_count").sum(),
                name=pl.lit("all"),
            )
            .select(
                "name",
                "a",
                "a_count",
                "b",
                "b_count",
                "diff",
                "diff_count",
            ),
        ]
    )
    df_cmp_sm_ab = (
        df_cmp_ab.group_by(["sample"])
        .agg(a=pl.col("length_a").sum(), b=pl.col("length_b").sum())
        .with_columns(diff=pl.col("b") - pl.col("a"))
    )

    # Break by AFR/Non-AFR
    draw_combined_err_diff(
        df_cmp_name_ab,
        output_prefix=f"{args.output_prefix}_all_count",
        name_colors=name_colors,
        by="count",
        include_all=True,
    )
    draw_combined_err_diff(
        df_cmp_name_ab,
        output_prefix=f"{args.output_prefix}_all_length",
        name_colors=name_colors,
        by="length",
        include_all=True,
    )
    draw_combined_err_diff(
        df_cmp_name_ab,
        output_prefix=f"{args.output_prefix}_count",
        name_colors=name_colors,
        by="count",
        include_all=False,
    )
    draw_combined_err_diff(
        df_cmp_name_ab,
        output_prefix=f"{args.output_prefix}_length",
        name_colors=name_colors,
        by="length",
        include_all=False,
    )
    median_diff_called_bases = df_cmp_sm_ab["diff"].median()
    perc_median_diff_called_bases = (
        median_diff_called_bases / DIPLOID_GENOME_SIZE
    ) * 100
    print(
        f"Median reduction in called bases (Mbp) from a to b: {median_diff_called_bases} ({perc_median_diff_called_bases}%)"
    )
    df_cmp_ab.write_csv(f"{args.output_prefix}_sm_name_diff.tsv", separator="\t")
    df_cmp_sm_ab.write_csv(f"{args.output_prefix}_sm_diff.tsv", separator="\t")
    df_cmp_name_ab.write_csv(f"{args.output_prefix}_name_diff.tsv", separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
