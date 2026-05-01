import os
import sys
import argparse
import matplotlib.colors

import polars as pl
import intervaltree as it
import matplotlib.pyplot as plt

from functools import lru_cache
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from collections import defaultdict
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm


BEDPE_IMPG_COLS = ("qchrom", "qst", "qend", "rchrom", "rst", "rend", "ritv")
TYPES = (
    "insertion",
    "deletion",
    "false_dup",
    "misjoin",
    "collapse",
    "mismatch",
    "scaffold",
    "other_repeat",
    "softclip",
)
CMAP = colormaps.get_cmap("OrRd")
NO_OVL_COLOR = "gray"
CTG_BRK_COLOR = "black"
NORM = matplotlib.colors.Normalize(vmin=0, vmax=100)
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)
SEGDUP_COLORS = dict(
    zip(
        [
            r"Less than 90% similarity",
            r"90 - 98% similarity",
            r"98 - 99% similarity",
            r"Greater than 99% similarity",
        ],
        ["#800080", "#808080", "#ffff00ff", "#ffa500"],
    )
)
plt.rcParams["font.family"] = "Arial"


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


def draw_perc_ovl_legend(fig: plt.Figure):
    # Draw perc overlap
    increment = 10
    bounds = [brkpt for brkpt in range(0, 110, increment)]

    # https://matplotlib.org/stable/users/explain/colors/colorbar_only.html#discrete-and-extended-colorbar-with-continuous-colorscale
    norm = BoundaryNorm(bounds, CMAP.N, extend="both")

    ax: Axes = fig.gca()
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=CMAP),
        cax=ax,
        orientation="horizontal",
    )

    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlabel("Percent NucFlag call overlap", fontsize=18)

    labels = ["No overlap", "Contig break"]
    handles = [
        Patch(edgecolor="black", facecolor=NO_OVL_COLOR),
        Patch(edgecolor="black", facecolor=CTG_BRK_COLOR),
    ]
    fig.legend(
        labels=labels,
        handles=handles,
        title_fontsize=18,
        fontsize=16,
        ncols=len(labels),
        loc="center",
        alignment="left",
        bbox_to_anchor=(0.5, -0.2),
        **LEGEND_KWARGS,
    )


@lru_cache
def convert_item_rgb(item_rgb: str) -> str:
    channels = [float(v) / 255.0 for v in item_rgb.split(",")]
    return matplotlib.colors.rgb2hex(channels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i", "--input_liftover", help="impg liftover BEDPE file", required=True
    )
    ap.add_argument(
        "-a", "--annotations", nargs="+", help="Annotations relative to reference."
    )
    ap.add_argument("-al", "--annotation_labels", nargs="+", help="Annotation labels.")
    ap.add_argument("-b", "--bed_calls", nargs="+", help="NucFlag bed files.")
    ap.add_argument("-l", "--labels", nargs="+", help="Labels for NucFlag calls.")
    ap.add_argument("-f", "--fais", nargs="+", help="Fasta indexes for each label.")
    ap.add_argument("-c", "--colors", nargs="+", help="Colors for NucFlag calls.")
    ap.add_argument("-o", "--output_dir", default=".", help="Output plot dir.")
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    label_colors: dict[str, str] = {}
    itrees_calls: defaultdict[str, it.IntervalTree] = defaultdict(it.IntervalTree)
    ctg_lengths: defaultdict[str, dict[str, int]] = defaultdict(dict)
    for label, color, calls, fais in zip(
        args.labels, args.colors, args.bed_calls, args.fais, strict=True
    ):
        df_calls = pl.read_csv(
            calls,
            separator="\t",
            columns=(0, 1, 2, 3),
            has_header=False,
            comment_prefix="#",
            new_columns=("chrom", "st", "end", "name"),
        ).filter(pl.col("name").is_in(TYPES))

        df_fais = pl.read_csv(
            fais,
            columns=[0, 1],
            separator="\t",
            has_header=False,
            new_columns=["chrom", "length"],
        )
        ctg_lengths[label] = dict(df_fais.iter_rows())
        label_colors[label] = color
        for row in df_calls.iter_rows(named=True):
            if row["st"] == row["end"]:
                continue
            itrees_calls[f"{label}_{row['chrom']}"].add(
                it.Interval(row["st"], row["end"], row["name"])
            )

    all_labels = set(label_colors.keys())
    annotations: defaultdict[str, defaultdict[str, set[it.Interval]]] = defaultdict(
        lambda: defaultdict(set)
    )
    dfs_annot = []
    for lbl, annot in zip(args.annotation_labels, args.annotations, strict=True):
        df_annot = pl.read_csv(
            annot,
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "st",
                "end",
                "name",
                "score",
                "strand",
                "tst",
                "tend",
                "item_rgb",
            ],
            columns=list(range(9)),
        ).with_columns(lbl=pl.lit(lbl))
        dfs_annot.append(df_annot)

        for row in df_annot.iter_rows(named=True):
            annotations[lbl][row["chrom"]].add(
                it.Interval(
                    row["st"],
                    row["end"],
                    (row["name"], convert_item_rgb(row["item_rgb"])),
                )
            )
    df_all_annot = pl.concat(dfs_annot)

    n_annotations = len(annotations)
    label_order = {lbl: i + n_annotations for i, lbl in enumerate(label_colors.keys())}

    df_liftover = (
        pl.read_csv(
            args.input_liftover,
            new_columns=list(BEDPE_IMPG_COLS),
            has_header=False,
            separator="\t",
        )
        .filter(pl.col("qchrom") != pl.col("rchrom"))
        .with_columns(
            rlabel=pl.col("rchrom").str.extract(r"^(.+)_[^_]*?$"),
            qlabel=pl.col("qchrom").str.extract(r"^(.+)_[^_]*?$"),
        )
        # No self-alignments
        .filter(pl.col("rlabel") != pl.col("qlabel"))
        .select(*BEDPE_IMPG_COLS, "rlabel", "qlabel")
    )
    figs = {}
    for rchrom in df_liftover["rchrom"].unique():
        fig, axes = plt.subplots(
            nrows=len(label_colors) + len(annotations),
            ncols=1,
            sharex=True,
            layout="constrained",
            figsize=(16, 4),
        )

        # Draw annotations
        for i, (lbl, chrom_annot) in enumerate(annotations.items()):
            annots = chrom_annot.get(rchrom)
            if not annots:
                continue
            ax_annot: Axes = axes[i]
            ax_annot.tick_params(axis="both", which="major", labelsize=16)
            minimalize_ax(ax_annot, remove_ticks=True)
            ax_annot.set_ylabel(lbl, rotation=0, ha="right", fontsize=18)
            for annot in annots:
                ax_annot.axvspan(
                    annot.begin, annot.end, color=annot.data[1], label=annot.data[0]
                )

        for lbl, idx in label_order.items():
            ax: Axes = axes[idx]
            ax.tick_params(axis="both", which="major", labelsize=16)
            lbl_color = label_colors[lbl]
            ax.set_ylabel(lbl, rotation=0, ha="right", color=lbl_color, fontsize=18)
            if idx == len(axes) - 1:
                ax.set_xlabel("Position (Mbp)", fontsize=18)
            ax.margins(x=0)
            ax.set_yticks([], [])
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

        figs[rchrom] = (fig, axes)

    for grp, df_grp in df_liftover.group_by(["ritv"]):
        rchrom, ritv = grp[0].split(":")
        rst, rend = ritv.split("-")
        rst, rend = int(rst), int(rend)
        print(f"On {rchrom}:{rst}-{rend}", file=sys.stderr)
        fig, axes = figs[rchrom]
        rlabel, _ = rchrom.rsplit("_")
        ax_ref: Axes = axes[label_order[rlabel]]

        itree_ref = itrees_calls.get(rchrom)
        if itree_ref:
            ovl_ref = itree_ref.overlap(rst, rend)
            if ovl_ref:
                ovl_length_ref = sum(
                    min(itv.end, rend) - max(itv.begin, rst) for itv in ovl_ref
                )
                prop_calls = (ovl_length_ref / (rend - rst)) * 100
                assert round(prop_calls) <= 100.0, (
                    f"Overlap intervals ({ovl_ref}, {ovl_length_ref}) exceeds ref region bounds. {prop_calls}"
                )

                color = CMAP(NORM(prop_calls))
                ax_ref.axvspan(rst, rend, color=color)

        # Color labels with no aligned region black.
        missing_labels = all_labels.difference([rlabel, *df_grp["qlabel"]])
        for label in missing_labels:
            i = label_order[label]
            ax: Axes = axes[i]
            ax.axvspan(rst, rend, color=NO_OVL_COLOR)

        merged_rows = []
        # Merge overlaps in query
        for _, df_grp_qchrom in df_grp.group_by(["qchrom"]):
            itree_grp_qchrom = it.IntervalTree(
                it.Interval(qrow["qst"], qrow["qend"], qrow)
                for qrow in df_grp_qchrom.iter_rows(named=True)
            )
            itree_grp_qchrom.merge_overlaps(
                data_reducer=lambda curr_data, new_data: (
                    {
                        "qchrom": curr_data["qchrom"],
                        "qst": curr_data["qst"],
                        "qend": curr_data["qend"],
                        "rchrom": curr_data["rchrom"],
                        "rst": min(curr_data["rst"], new_data["rst"]),
                        "rend": max(curr_data["rend"], new_data["rend"]),
                        "ritv": curr_data["ritv"],
                        "rlabel": curr_data["rlabel"],
                        "qlabel": curr_data["qlabel"],
                    }
                )
            )
            merged_rows.extend(
                (
                    {**mdata, "qst": st, "qend": end}
                    for st, end, mdata in itree_grp_qchrom.iter()
                )
            )

        df_grp_merged = pl.DataFrame(merged_rows, orient="row")

        for row in df_grp_merged.iter_rows(named=True):
            itree = itrees_calls.get(row["qchrom"])
            if not itree:
                continue

            ovl = itree.overlap(row["qst"], row["qend"])
            if not ovl:
                continue

            i = label_order[row["qlabel"]]
            ax: Axes = axes[i]
            window_length = row["rend"] - row["rst"]
            # Clip to region
            ovl_length = sum(
                min(itv.end, row["qend"]) - max(itv.begin, row["qst"]) for itv in ovl
            )
            # Clamp to 100
            prop_calls = min(100.0, (ovl_length / window_length) * 100)
            # Mark if at contig end
            ctg_end = ctg_lengths[row["qlabel"]][row["qchrom"]]
            at_ctg_end = any(itv.begin == 0 or itv.end == ctg_end for itv in ovl)
            if at_ctg_end:
                color = CTG_BRK_COLOR
            else:
                color = CMAP(NORM(prop_calls))
            ax.axvspan(row["rst"], row["rend"], color=color)

    for rchrom, (fig, axes) in figs.items():
        for ax in axes:
            ax.xaxis.set_major_formatter(formatter=lambda x, _: f"{x / 1_000_000:.1f}")

        fig.suptitle(
            rchrom.replace("CHM13v2.0_chr", "Chromosome "),
            x=0.0,
            fontsize=24,
            weight="bold",
            horizontalalignment="left",
            verticalalignment="center",
        )

        outfile_prefix = os.path.join(args.output_dir, f"{rchrom}")
        fig.savefig(f"{outfile_prefix}.png", bbox_inches="tight", dpi=300)
        plt.close(fig)

        # Create figure for each
        fig_legend, axes_legend = plt.subplots(
            nrows=len(args.annotation_labels), layout="constrained"
        )

        for i, lbl in enumerate(args.annotation_labels):
            ax_legend: Axes = axes_legend[i]
            minimalize_ax(ax_legend, remove_ticks=True)
            if lbl == "Segmental duplications":
                # Use predefined colors. May include those not shown in data.
                ax_legend.legend(
                    labels=SEGDUP_COLORS.keys(),
                    handles=[
                        Patch(edgecolor="black", facecolor=color)
                        for color in SEGDUP_COLORS.values()
                    ],
                    loc="center",
                    title=lbl,
                    title_fontsize=18,
                    fontsize=16,
                    ncols=len(SEGDUP_COLORS.keys()),
                    alignment="left",
                    **LEGEND_KWARGS,
                )
            else:
                df_annot = df_all_annot.filter(
                    pl.col("chrom").eq(pl.lit(rchrom)) & pl.col("lbl").eq(lbl)
                )
                if lbl == "Satellite structure":
                    # Only take first part of censat name dhor_*
                    uniq_elems = {
                        name.split("_")[0]: convert_item_rgb(item_rgb)
                        for name, item_rgb in df_annot.select("name", "item_rgb")
                        .unique()
                        .sort(by="name")
                        .iter_rows()
                    }
                else:
                    uniq_elems = {
                        name: convert_item_rgb(item_rgb)
                        for name, item_rgb in df_annot.select("name", "item_rgb")
                        .unique()
                        .sort(by="name")
                        .iter_rows()
                    }
                ax_legend.legend(
                    title=lbl,
                    labels=uniq_elems.keys(),
                    loc="center",
                    handles=[
                        Patch(edgecolor="black", facecolor=color)
                        for color in uniq_elems.values()
                    ],
                    title_fontsize=18,
                    fontsize=16,
                    ncols=len(uniq_elems),
                    alignment="left",
                    **LEGEND_KWARGS,
                )

        fig_legend.savefig(f"{outfile_prefix}_legend.png", bbox_inches="tight", dpi=300)
        fig_legend.savefig(f"{outfile_prefix}_legend.pdf", bbox_inches="tight", dpi=300)
        plt.close(fig_legend)

    fig_ovl_legend = plt.figure(figsize=(8, 1), layout="constrained")
    draw_perc_ovl_legend(fig_ovl_legend)
    output_prefix_ovl = os.path.join(args.output_dir, "ovl_legend")
    fig_ovl_legend.savefig(f"{output_prefix_ovl}.png", bbox_inches="tight", dpi=300)
    fig_ovl_legend.savefig(f"{output_prefix_ovl}.pdf", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
