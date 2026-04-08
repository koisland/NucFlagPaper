import os
import sys
import argparse
import matplotlib.colors

import polars as pl
import intervaltree as it
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from collections import defaultdict
from functools import lru_cache


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
CMAP = plt.get_cmap("OrRd")
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


def draw_perc_ovl_legend(obj: Axes | plt.Figure):
    # Draw perc overlap
    legend_handles_perc_calls = {}
    increment = 10
    for brkpt in range(0, 110, increment):
        color = CMAP(NORM(brkpt))
        legend_handles_perc_calls[f"{brkpt}%"] = Patch(
            edgecolor="black", facecolor=color
        )

    obj.legend(
        labels=[*legend_handles_perc_calls.keys(), "No overlap", "Contig break"],
        handles=[
            *legend_handles_perc_calls.values(),
            Patch(edgecolor="black", facecolor=NO_OVL_COLOR),
            Patch(edgecolor="black", facecolor=CTG_BRK_COLOR),
        ],
        title="Percent NucFlag call overlap",
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
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
            itrees_calls[f"{label}_{row['chrom']}"].add(
                it.Interval(row["st"], row["end"], row["name"])
            )

    all_labels = set(label_colors.keys())
    annotations: defaultdict[str, defaultdict[str, set[it.Interval]]] = defaultdict(
        lambda: defaultdict(set)
    )
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
        )

        for row in df_annot.iter_rows(named=True):
            annotations[lbl][row["chrom"]].add(
                it.Interval(
                    row["st"],
                    row["end"],
                    (row["name"], convert_item_rgb(row["item_rgb"])),
                )
            )

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
            figsize=(16, 6),
        )

        # Draw annotations
        for i, (lbl, chrom_annot) in enumerate(annotations.items()):
            annots = chrom_annot.get(rchrom)
            if not annots:
                continue
            ax_annot: Axes = axes[i]
            minimalize_ax(ax_annot, remove_ticks=True)
            ax_annot.set_ylabel(lbl, rotation=0, ha="right")
            for annot in annots:
                ax_annot.axvspan(
                    annot.begin, annot.end, color=annot.data[1], label=annot.data[0]
                )

        for lbl, idx in label_order.items():
            ax: Axes = axes[idx]
            lbl_color = label_colors[lbl]
            ax.set_ylabel(lbl, rotation=0, ha="right", color=lbl_color)
            if idx == len(axes) - 1:
                ax.set_xlabel("Position (Mbp)")
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

        draw_perc_ovl_legend(fig)
        fig.savefig(
            os.path.join(args.output_dir, f"{rchrom}.png"), bbox_inches="tight", dpi=300
        )
        plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
