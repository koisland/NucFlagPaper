import sys
import random
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from dataclasses import dataclass
from collections import Counter, defaultdict
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib import ticker
from intervaltree import Interval, IntervalTree

plt.rcParams["font.family"] = "Arial"

LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    loc="center left",
    alignment="left",
)
ANNOT_COLS = {
    "column_9": "qstrand",
    "column_11": "qchrom",
    "column_12": "qst",
    "column_13": "qend",
    "column_14": "name",
    "column_16": "strand",
    "column_19": "item_rgb",
    "column_4": "rchrom",
    "column_5": "rst",
    "column_6": "rend",
}
DUPMASKER_COLORS = {
    "chr1": "#ffe4e1",
    "chr2": "#000000",
    "chr3": "#2f4f4f",
    "chr4": "#191970",
    "chr5": "#6495ed",
    "chr6": "#0000cd",
    "chr7": "#00bfff",
    "chr8": "#4682b4",
    "chr9": "#00ffff",
    "chr10": "#e0ffff",
    "chr11": "#013220",
    "chr12": "#00ff7f",
    "chr13": "#7cfc00",
    "chr14": "#ffff00",
    "chr15": "#ffd700",
    "chr16": "#bc8f8f",
    "chr17": "#cd5c5c",
    "chr18": "#a52a2a",
    "chr19": "#ff6347",
    "chr20": "#cd853f",
    "chr21": "#ff0000",
    "chr22": "#ff69b4",
    "chrX": "#ff1493",
    "chrY": "#800000",
    "chr_random": "#9400d3",
}
# https://humanpangenome.org/samples/
AFR_POP_LABELS = {"ACB", "GWD", "ESN", "MKK", "LWK", "ASL", "YRI", "MSL", "ASW"}
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
STRUCTURAL_ERRORS = {"insertion", "deletion", "other_repeat", *LARGE_ERRORS}
RELEASE_COLORS = {"Release 1": "red", "Release 2": "blue"}
NON_ERROR_CALLS = {"correct", "het_or_mismap"}


@dataclass
class OverlapCounts:
    errors: int = 0
    breaks: int = 0


def item_rgb_to_hex(color: str) -> str:
    return to_hex(tuple(int(c) / 255.0 for c in color.split(",")))


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False) -> None:
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
            labelbottom=False,
        )


def make_reorient_relative_df(df: pl.DataFrame) -> pl.DataFrame:
    """
    Reorient intervals based on CHM13 alignment.
    Also make relative. Requires all annotation data to be present as takes global minimum start position.
    """
    df_out = (
        df.with_columns(
            qst=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qend"))
            .otherwise(pl.col("qst")),
            qend=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qst"))
            .otherwise(pl.col("qend")),
            strand=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(
                pl.when(pl.col("strand").eq(pl.lit("-")))
                .then(pl.lit("+"))
                .when(pl.col("strand").eq(pl.lit("+")))
                .then(pl.lit("-"))
                .otherwise(pl.col("strand"))
            )
            .otherwise(pl.col("strand")),
        )
        .with_columns(
            st=(pl.col("qst") - pl.col("qst").min().over("qchrom")) + pl.col("rst"),
            end=(pl.col("qend") - pl.col("qst").min().over("qchrom")) + pl.col("rst"),
            mtch=pl.col("qchrom").str.extract_groups(r"^(?<sm>.*?)#(?<hap>1|2)#"),
        )
        .with_columns(
            pl.col("st") - pl.col("st").min(),
            pl.col("end") - pl.col("st").min(),
        )
        .unnest("mtch")
        .sort(["rchrom", "hap", "st"])
    )
    # # Debugging
    # df_out_sm = df_out.filter(pl.col("qchrom").str.contains("HG01243#2"))
    # fig, ax = plt.subplots()
    # ax: Axes

    # for row in df_out_sm.iter_rows(named=True):
    #     ax.axvspan(row["st"], row["end"], color=item_rgb_to_hex(row["item_rgb"]))
    # fig.savefig("test.png")
    # df_out_sm.write_csv("test.csv")
    # breakpoint()
    return df_out


def draw_r1_r2_smn(
    df_all: pl.DataFrame,
    nucflag_colors: dict[str, str],
    output_prefix: str,
    figsize: tuple[int, int] = (16, 24),
):
    sm_haps = sorted(set(df_all.select("sm", "hap").iter_rows()))
    nrows = (len(sm_haps) * 2) + 1
    ht_ratios = [0.25, 1.0] * len(sm_haps)
    ht_ratios.append(0.25)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=2,
        layout="constrained",
        sharex=True,
        height_ratios=ht_ratios,
        figsize=figsize,
    )

    for col in (0, 1):
        ax_bottom: Axes = axes[-1, col]
        for spine in ["top", "left", "right"]:
            ax_bottom.spines[spine].set_visible(False)
        ax_bottom.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
        )
        ax_bottom.xaxis.set_major_formatter(lambda x, pos: f"{x / 1_000_000:.1f}")
        ax_bottom.set_xlabel("Relative position (Mbp)")

    breaks_counter = defaultdict(Counter)
    breaks_color = {True: "#000000", False: "#72bcd4"}
    breaks_labels = {True: "Break", False: "Complete"}
    # Mark if has error.
    df_all = df_all.with_columns(
        has_error=pl.col("name")
        .is_in(STRUCTURAL_ERRORS)
        .any()
        .over(["sm", "hap", "release"])
    )
    err_norm_itree: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    break_norm_itree: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    median_length_release: dict[str, float] = dict(
        df_all.filter(pl.col("dtype").eq(pl.lit("dupmasker")))
        .group_by(["release"])
        .agg((pl.col("end").max() - pl.col("st").min()).median())
        .iter_rows()
    )
    for row_offset, annot in enumerate(("nucflag", "dupmasker")):
        df_annot = df_all.filter(pl.col("dtype").eq(pl.lit(annot)))
        for row, (sm, hap) in enumerate(sm_haps):
            row = (row * 2) + row_offset

            for col, release in enumerate(("Release 1", "Release 2")):
                median_length = median_length_release[release]
                df_sm_hap_annot = df_annot.filter(
                    pl.col("sm").eq(pl.lit(sm))
                    & pl.col("hap").eq(pl.lit(hap))
                    & pl.col("release").eq(pl.lit(release))
                )
                # max_end_position = df_sm_hap_annot["end"].max()
                ax: Axes = axes[row, col]
                ax.set_ylim(0, 1)
                print(row, col, file=sys.stderr)

                ax.margins(x=0, y=0)

                minimalize_ax(ax, remove_ticks=True)

                polygons = []
                ht = ax.get_ylim()[1]
                for itv in df_sm_hap_annot.iter_rows(named=True):
                    hexcolor = item_rgb_to_hex(itv["item_rgb"])
                    if annot == "dupmasker":
                        # Draw triangle
                        # https://matplotlib.org/devdocs/api/_as_gen/matplotlib.patches.Polygon.html
                        if itv["strand"] == "+":
                            vertices = [
                                (itv["st"], 0),
                                (itv["st"], ht),
                                (itv["end"], ht / 2),
                            ]
                        else:
                            vertices = [
                                (itv["end"], 0),
                                (itv["end"], ht),
                                (itv["st"], ht / 2),
                            ]
                        poly = Polygon(xy=vertices, color=hexcolor)
                        polygons.append(poly)
                    else:
                        if itv["name"] in NON_ERROR_CALLS:
                            continue
                        err_norm_itree[release].add(Interval(itv["st"], itv["end"]))

                        # color = hexcolor if itv["name"] in STRUCTURAL_ERRORS else "#000000"
                        color = "#000000"
                        ax.axvspan(itv["st"], itv["end"], color=color)

                if annot == "dupmasker":
                    # # Also draw boundary rect
                    # # NOTE: This is not good for visualization because it obfuscates misaligned regions and gaps.
                    # min_st, max_end = df_sm_hap_annot["st"].min(), df_sm_hap_annot["end"].max()
                    # ax.axvspan(min_st, max_end, color="#808080", alpha=0.5)
                    ax.add_collection(PatchCollection(polygons, match_original=True))

                    has_break_or_error = df_sm_hap_annot["has_error"].first() is True
                    for _, df_grp in df_sm_hap_annot.group_by(["qchrom"]):
                        # Is a break since if close to end of contig. Draw dashed line.
                        end_len_diff = abs(
                            df_grp["qend"].max() - df_grp["length"].first()
                        )
                        if df_grp["qst"].min() < 50000:
                            x_line = df_grp["st"].min()
                        elif end_len_diff < 50000:
                            x_line = df_grp["end"].max()
                        else:
                            continue
                        print(row, col, "break")
                        # norm_xst = (x_line / max_end_position) * median_length
                        break_norm_itree[release].add(Interval(x_line, x_line + 1))
                        ax.axvline(
                            x=x_line,
                            ymin=0,
                            ymax=ht,
                            linestyle="dashed",
                            linewidth=1,
                            color="black",
                            zorder=2,
                        )
                        has_break_or_error = True

                    breaks_counter[release][has_break_or_error] += 1
                    color = breaks_color[has_break_or_error]
                    ax.set_ylabel(
                        f"{sm}_hap{hap}",
                        color=color,
                        rotation=0,
                        ha="right",
                        va="center",
                        fontsize=12,
                    )

    # Draw relative position of breaks
    fig_cov, axes_cov = plt.subplots(
        ncols=2,
        nrows=1,
        layout="constrained",
        figsize=(figsize[0], 2),
        sharex=True,
        sharey=True,
    )
    window = 5000
    for i, release in enumerate(("Release 1", "Release 2")):
        ax: Axes = axes_cov[i]
        ax.margins(x=0)
        breaks_itree = break_norm_itree[release]
        errors_itree = err_norm_itree[release]
        median_length = median_length_release[release]

        num, rem = divmod(median_length, window)
        num = int(num)
        final_start = num * window
        # All intervals of window size
        itvs = [
            Interval((i - 1) * window, i * window, OverlapCounts())
            for i in range(1, num + 1)
        ]
        itvs.append(Interval(final_start, final_start + rem, OverlapCounts()))

        for itv in itvs:
            cnts: OverlapCounts = itv.data
            for bitv in breaks_itree.overlap(itv):
                cnts.breaks += 1

            for eitv in errors_itree.overlap(itv):
                cnts.errors += 1

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        ax.xaxis.set_major_formatter(lambda x, pos: f"{x / 1_000_000:.1f}")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(500_000))
        ax.set_xlabel("Relative position (Mbp)")
        ax.set_ylabel("Cumulative # of calls")
        ax.set_title(release, color=RELEASE_COLORS[release])

        # Then draw the bars
        xpos, breaks, errors = [], [], []
        for itv in itvs:
            xpos.append(itv.begin)
            errors.append(itv.data.errors)
            breaks.append(itv.data.breaks)

        ax.bar(x=xpos, height=errors, width=window, color="black", label="Errors")
        ax.bar(x=xpos, height=breaks, width=window, color="orange", label="Breaks")
        ax.legend(**LEGEND_KWARGS | {"loc": "upper right"})

    fig_cov.savefig(
        f"{output_prefix}_relpos_breaks_errors.png", bbox_inches="tight", dpi=300
    )
    fig_cov.savefig(
        f"{output_prefix}_relpos_breaks_errors.pdf", bbox_inches="tight", dpi=300
    )

    # Draw breakdown of errors in locus
    fig_nucflag, ax_nucflag = plt.subplots(layout="constrained", figsize=(12, 4))
    ax_nucflag: Axes
    df_nucflag = (
        df_all.filter(pl.col("dtype").eq(pl.lit("nucflag")))
        .group_by(["release", "name"])
        .agg(length=(pl.col("end") - pl.col("st")).sum(), count=pl.col("end").count())
        .sort(by=["name", "release"])
    )
    sns.barplot(
        data=df_nucflag,
        x="name",
        y="count",
        hue="release",
        order=[*LARGE_ERRORS, *SMALL_ERRORS],
        hue_order=["Release 1", "Release 2"],
        palette=RELEASE_COLORS,
        ax=ax_nucflag,
    )
    for cont in ax_nucflag.containers:
        ax_nucflag.bar_label(
            cont,
            fontsize=12,
            label_type="edge",
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        )
    ax_nucflag.tick_params(axis="both", which="major", labelsize=14)
    # Draw median
    df_median_count = df_nucflag.group_by(["release"]).agg(pl.col("count").median())
    release_to_median = {
        release: median for release, median in df_median_count.iter_rows()
    }
    ax_nucflag_yticks, ax_nucflag_yticklabels = (
        list(ax_nucflag.get_yticks()),
        ax_nucflag.get_yticklabels(),
    )

    # Map median to release
    median_to_release = {}
    for release, median in release_to_median.items():
        ax_nucflag.axhline(y=median, color=RELEASE_COLORS[release], linestyle="dotted")
        ax_nucflag_yticks.append(median)
        ax_nucflag_yticklabels.append(f"{median:.0f}")
        median_to_release[f"{median:.0f}"] = release
    ax_nucflag.set_yticks(ax_nucflag_yticks, ax_nucflag_yticklabels)

    # Then color them
    for lbl in ax_nucflag.get_yticklabels():
        release = median_to_release.get(f"{float(lbl.get_text()):.0f}")
        if not release:
            continue
        color = RELEASE_COLORS[release]
        lbl.set_color(color)
        lbl.set_path_effects([pe.withStroke(linewidth=1, foreground="white")])

    ax_nucflag.set_xlabel(None)
    ax_nucflag.set_ylabel("# of calls", fontsize=14)
    for lbl in ax_nucflag.xaxis.get_majorticklabels():
        lbl.set_rotation(45)
        lbl.set_color(nucflag_colors[lbl.get_text()])
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")
        lbl.set_path_effects([pe.withStroke(linewidth=0.5, foreground="black")])

    for spine in ("top", "right"):
        ax_nucflag.spines[spine].set_visible(False)

    sns.move_legend(
        ax_nucflag, fontsize=14, title=None, **LEGEND_KWARGS | {"loc": "upper left"}
    )

    fig_nucflag.savefig(f"{output_prefix}_errors.pdf", bbox_inches="tight", dpi=300)
    fig_nucflag.savefig(f"{output_prefix}_errors.png", bbox_inches="tight", dpi=300)

    # Draw hbar of number with breaks
    fig_breaks, axes_breaks = plt.subplots(
        layout="constrained", ncols=2, nrows=1, figsize=(figsize[0], 0.75), sharex=True
    )
    for i, ax_breaks in enumerate(axes_breaks):
        ax_breaks: Axes
        release = f"Release {i + 1}"
        sns.barplot(
            x=[sum(breaks_counter[release].values())],
            y=[release],
            color=breaks_color[False],
            ax=ax_breaks,
            label=breaks_labels[False],
            legend=None,
        )
        sns.barplot(
            x=[breaks_counter[release][True]],
            y=[release],
            color=breaks_color[True],
            ax=ax_breaks,
            label=breaks_labels[True],
            legend=None,
        )
        # NOTE: This will be centered so needs to be adjusted in post
        for cont in ax_breaks.containers:
            ax_breaks.bar_label(
                cont,
                label_type="center",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
            )

        for label in ax_breaks.get_yticklabels():
            label.set_color(RELEASE_COLORS[release])

        for spine in ("top", "left", "right"):
            ax_breaks.spines[spine].set_visible(False)

        ax_breaks.set_xlabel("# complete and structurally error-free")

        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
        )

    fig_breaks.savefig(f"{output_prefix}_breaks.pdf", bbox_inches="tight", dpi=300)
    fig_breaks.savefig(f"{output_prefix}_breaks.png", bbox_inches="tight", dpi=300)

    # NOTE: Needs to be manually moved to same x as ylbl
    ax_r1_top: Axes = axes[0, 0]
    ax_r1_top.set_title(
        "SMN1/2 locus (Release 1)",
        color=RELEASE_COLORS["Release 1"],
        loc="left",
        fontsize=18,
    )
    ax_r2_top: Axes = axes[0, 1]
    ax_r2_top.set_title(
        "SMN1/2 locus (Release 2)",
        color=RELEASE_COLORS["Release 2"],
        loc="left",
        fontsize=18,
    )

    fig.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1_dupmasker", nargs="+", type=str)
    ap.add_argument("--r2_dupmasker", nargs="+", type=str)
    ap.add_argument("--r1_nucflag", nargs="+", type=str)
    ap.add_argument("--r2_nucflag", nargs="+", type=str)
    ap.add_argument("--r1_fai", nargs="+", type=str)
    ap.add_argument("--r2_fai", nargs="+", type=str)
    ap.add_argument("--n_subset", default=10, type=int)
    ap.add_argument("--seed", default=7, type=int)
    ap.add_argument("--output_prefix", default="out", type=str)
    args = ap.parse_args()

    dfs_r1_nucflag = []
    dfs_r2_nucflag = []
    dfs_r1_dupmasker = []
    dfs_r2_dupmasker = []

    # contig lengths
    all_lengths = []
    for lbl, fais in (("Release 1", args.r1_fai), ("Release 2", args.r2_fai)):
        for file in fais:
            df_fai = pl.read_csv(
                file,
                separator="\t",
                columns=[0, 1],
                new_columns=["chrom", "length"],
                has_header=False,
            )
            for chrom, length in df_fai.iter_rows():
                all_lengths.append([lbl, chrom, length])

    df_all_lengths = pl.DataFrame(
        all_lengths, schema=["release", "qchrom", "length"], orient="row"
    )

    for file in args.r1_nucflag:
        df_r1_nucflag = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
            )
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r1_nucflag.append(df_r1_nucflag)

    for file in args.r2_nucflag:
        df_r2_nucflag = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
            )
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r2_nucflag.append(df_r2_nucflag)

    for file in args.r1_dupmasker:
        df_r1_dupmasker = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r1_dupmasker.append(df_r1_dupmasker)

    for file in args.r2_dupmasker:
        df_r2_dupmasker = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r2_dupmasker.append(df_r2_dupmasker)

    # Add length
    df_r1_nucflag = (
        pl.concat(dfs_r1_nucflag)
        .filter(pl.col("name").ne(pl.lit("correct")))
        .join(
            df_all_lengths.filter(pl.col("release").eq(pl.lit("Release 1"))),
            on="qchrom",
            how="left",
        )
    )
    df_r2_nucflag = (
        pl.concat(dfs_r2_nucflag)
        .filter(pl.col("name").ne(pl.lit("correct")))
        .join(
            df_all_lengths.filter(pl.col("release").eq(pl.lit("Release 1"))),
            on="qchrom",
            how="left",
        )
    )
    df_r1_dupmasker = pl.concat(dfs_r1_dupmasker).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("Release 1"))),
        on="qchrom",
        how="left",
    )
    df_r2_dupmasker = pl.concat(dfs_r2_dupmasker).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("Release 2"))),
        on="qchrom",
        how="left",
    )

    # Take min and max of same interval in case of alignment gap.
    # Also reorient if not forward oriented
    df_all = pl.concat(
        [
            df_r1_nucflag.with_columns(
                release=pl.lit("Release 1"), dtype=pl.lit("nucflag")
            ),
            df_r2_nucflag.with_columns(
                release=pl.lit("Release 2"), dtype=pl.lit("nucflag")
            ),
            df_r1_dupmasker.with_columns(
                release=pl.lit("Release 1"), dtype=pl.lit("dupmasker")
            ),
            df_r2_dupmasker.with_columns(
                release=pl.lit("Release 2"), dtype=pl.lit("dupmasker")
            ),
        ]
    )
    df_all = make_reorient_relative_df(df_all)

    nucflag_colors = {
        name: item_rgb_to_hex(color)
        for name, color in df_all.filter(pl.col("dtype").eq("nucflag"))
        .select("name", "item_rgb")
        .unique()
        .iter_rows()
    }

    # Draw subset of 10
    random.seed(args.seed)
    subset_sm = random.sample(list(df_all["sm"].unique()), args.n_subset)
    df_subset_all = df_all.filter(pl.col("sm").is_in(subset_sm))

    draw_r1_r2_smn(
        df_all=df_all, nucflag_colors=nucflag_colors, output_prefix=args.output_prefix
    )
    draw_r1_r2_smn(
        df_all=df_subset_all,
        nucflag_colors=nucflag_colors,
        output_prefix=f"{args.output_prefix}_subset{args.n_subset}",
        figsize=(16, 8),
    )

    fig_legend, axes_legend = plt.subplots(nrows=3, ncols=1, layout="tight")
    axes_legend: list[Axes]
    for ax in axes_legend:
        minimalize_ax(ax, remove_ticks=True)
    axes_legend[0].legend(
        handles=[
            Patch(facecolor=color, edgecolor="black")
            for lbl, color in nucflag_colors.items()
        ],
        labels=nucflag_colors.keys(),
        ncols=len(nucflag_colors) // 2,
        title="NucFlag calls",
        **LEGEND_KWARGS,
    )
    axes_legend[1].legend(
        handles=[
            Patch(facecolor=color, edgecolor="black")
            for lbl, color in DUPMASKER_COLORS.items()
        ],
        labels=DUPMASKER_COLORS.keys(),
        ncols=len(DUPMASKER_COLORS) // 2,
        title="Segmental duplications",
        **LEGEND_KWARGS,
    )
    # https://www.geeksforgeeks.org/python/custom-legends-with-matplotlib/
    axes_legend[2].legend(
        handles=[Line2D([0], [0], color="black", linestyle="dashed")],
        labels=["Contig break"],
        **LEGEND_KWARGS | {"handlelength": 2.0},
    )
    fig_legend.savefig(f"{args.output_prefix}_legend.png", bbox_inches="tight", dpi=300)
    fig_legend.savefig(f"{args.output_prefix}_legend.pdf", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
