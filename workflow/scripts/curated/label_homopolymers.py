import os
import sys
import pyfaidx
import argparse

import numpy as np
import polars as pl
import intervaltree as it
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from collections import defaultdict, Counter
from concurrent.futures import Future, ProcessPoolExecutor, as_completed

NT_COLORS = {
    "A": "#FF0000",
    "T": "#0000FF",
    "G": "#009600",
    "C": "#D17105",
    # False positive
    "N": "#000000",
}
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)


def minimalize_ax(
    ax: Axes,
    *,
    spines: tuple[str, ...] = ("left", "right", "bottom", "top"),
    remove_ticks: bool = False,
) -> None:
    for spine in spines:
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


def draw_stacked_bars(
    ax: Axes,
    hp_bins: defaultdict[str, Counter[int, int]],
    override_color: str | None = None,
    override_legend_label: str | None = None,
):
    max_homopolymer_length = max(max(lns.keys()) for _, lns in hp_bins.items())
    all_nt_homopolymer_cnts = np.zeros(max_homopolymer_length + 1)
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/bar_stacked.html
    for nt, cnt in hp_bins.items():
        if override_color:
            color = override_color
        else:
            color = NT_COLORS[nt]

        nt_homopolymer_ln = np.arange(max_homopolymer_length + 1)
        nt_homopolymer_cnts = np.zeros(max_homopolymer_length + 1)
        for h_ln, h_cnt in cnt.items():
            # If non-zero, is undefined bin.
            h_ln = max(h_ln - 1, 0)
            nt_homopolymer_cnts[h_ln] = h_cnt

        if override_legend_label:
            label = override_legend_label
        else:
            label = nt

        ax.bar(
            x=nt_homopolymer_ln,
            height=nt_homopolymer_cnts,
            color=color,
            alpha=0.5,
            label=label,
            bottom=all_nt_homopolymer_cnts,
        )
        all_nt_homopolymer_cnts += nt_homopolymer_cnts


def plot_homopolymers_stacked(
    homopolymer_bins: defaultdict[str, Counter[int, int]],
    outfile: str,
    xlabel: str = "Homopolymer length (bp)",
    ylabel: str = "Number of NucFlag homopolymer calls",
    other_homopolymer_bins: defaultdict[str, Counter[int, int]] | None = None,
    other_homopolymer_bin_color: str | None = None,
    other_homopolymer_bin_label: str | None = None,
):
    fig, ax = plt.subplots(layout="constrained", figsize=(8, 8))
    ax: Axes
    draw_stacked_bars(ax, homopolymer_bins)

    if other_homopolymer_bins:
        draw_stacked_bars(
            ax,
            other_homopolymer_bins,
            override_color=other_homopolymer_bin_color,
            override_legend_label=other_homopolymer_bin_label,
        )

    minimalize_ax(ax, spines=("top", "right"))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    handles, labels = ax.get_legend_handles_labels()
    labels_handles = {}
    for label, handle in zip(labels, handles):
        labels_handles[label] = handle

    ax.legend(
        labels=labels_handles.keys(),
        handles=labels_handles.values(),
        loc="upper right",
        **LEGEND_KWARGS,
    )
    fig.savefig(outfile, bbox_inches="tight", dpi=300)


def plot_homopolymers_split(
    homopolymer_bins: defaultdict[str, Counter[int, int]],
    outfile: str,
    xlabel: str = "Homopolymer length (bp)",
    ylabel: str = "Number of NucFlag homopolymer calls",
    *,
    include_n: bool,
):
    fig, axes = plt.subplots(
        layout="constrained",
        figsize=(16, 4),
        nrows=1,
        ncols=len(NT_COLORS.keys()) - int(not include_n),
        sharey=True,
    )
    # Draw all homopolymers as hatched bar
    max_homopolymer_length = max(max(lns.keys()) for _, lns in homopolymer_bins.items())

    for i, (nt, cnt) in enumerate(homopolymer_bins.items()):
        color = NT_COLORS[nt]
        ax: Axes = axes[i]
        nt_homopolymer_ln = np.arange(max_homopolymer_length + 1)
        nt_homopolymer_cnts = np.zeros(max_homopolymer_length + 1)
        for h_ln, h_cnt in cnt.items():
            # If non-zero, is undefined bin.
            h_ln = max(h_ln - 1, 0)
            nt_homopolymer_cnts[h_ln] = h_cnt

        # color 0 bin black
        colors = ["black"]
        colors.extend(color for _ in range(max_homopolymer_length))
        ax.bar(x=nt_homopolymer_ln, height=nt_homopolymer_cnts, color=colors, alpha=0.5)

        minimalize_ax(ax, spines=("top", "right"))
        ax.set_title(nt)

    labels_handles = {
        nt if nt != "N" else "Undefined": Patch(color=NT_COLORS[nt], alpha=0.5)
        for nt in ["N", *homopolymer_bins.keys()]
    }
    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    # https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
    fig.legend(
        labels=labels_handles.keys(),
        handles=labels_handles.values(),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(labels_handles.keys()),
        **LEGEND_KWARGS,
    )
    fig.savefig(outfile, bbox_inches="tight", dpi=300)


def process_one(fa: str, df_lcr: pl.DataFrame) -> tuple[str, it.IntervalTree]:
    fh_fai = pyfaidx.Faidx(fa)
    itree = it.IntervalTree()
    chrom = str(df_lcr["chrom"].first())
    print(f"On {chrom}...", file=sys.stderr)
    for chrom, st, end in df_lcr.iter_rows():
        # 0-based to 1-based (pyfaidx) for fetching.
        st, end = int(st) + 1, int(end)
        seq = fh_fai.fetch(chrom, st, end)
        nts = set(str(seq))
        # print(chrom, st, end, seq)
        is_homopolymer = len(nts) == 1
        if not is_homopolymer:
            continue
        nt = next(iter(nts))
        # Revert back to longdust coordinates 0-based.
        itree.add(it.Interval(st - 1, end, nt))
    print(f"Done with {chrom}...", file=sys.stderr)
    return chrom, itree


def main():  #
    ap = argparse.ArgumentParser()
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument(
        "-i",
        "--input_lcr_bed",
        type=argparse.FileType("rb"),
        help="Input low-complexity regions BED3 from longdust.",
    )
    grp.add_argument(
        "--homopolymers",
        type=argparse.FileType("rb"),
        help="Homopolymer BED4 file. Mutually exclusive with --input_lcr_bed With nt as name column",
    )
    ap.add_argument(
        "-r",
        "--reference",
        help="Reference fasta which was used with longdust.",
    )
    ap.add_argument(
        "-c",
        "--calls",
        type=argparse.FileType("rt"),
        required=True,
        help="NucFlag calls.",
    )
    ap.add_argument("-o", "--outfile", default=sys.stdout, help="All homopolymers.")
    ap.add_argument(
        "-x",
        "--outfile_incorrect",
        default=None,
        type=argparse.FileType("wt"),
        help="All homopolymers from NucFlag that have no overlap with longdust homopolymers.",
    )
    ap.add_argument(
        "-p",
        "--plot_prefix",
        default="homopolymers",
        help="Plot prefix for plots.",
    )

    args = ap.parse_args()

    df_nucflag_homopolymer_calls = pl.read_csv(
        args.calls, separator="\t", has_header=True
    ).filter(pl.col("name") == "homopolymer")
    all_homopolymer_bins = defaultdict(Counter)

    if not args.homopolymers:
        reference: str | None = args.reference
        if not reference:
            raise ValueError("Need reference.")

        df_bed_lcr = pl.read_csv(
            args.input_lcr_bed,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "st", "end"],
        )
        dfs_bed_lcr = df_bed_lcr.partition_by(["chrom"])

        homopolymers: defaultdict[str, it.IntervalTree] = defaultdict()
        with ProcessPoolExecutor(max_tasks_per_child=1, max_workers=12) as pool:
            futures: list[Future] = []
            for df in dfs_bed_lcr:
                futures.append(pool.submit(process_one, reference, df))

            for future in as_completed(futures):
                try:
                    chrom, df = future.result()
                    homopolymers[chrom] = df
                except Exception as res:
                    print(f"Failed for reason: {res}", file=sys.stderr)
                    continue

        print("Homopolymers filtered from longdust bed.", file=sys.stderr)
        itv: it.Interval
        for chrom, itree in homopolymers.items():
            for itv in itree.iter():
                all_homopolymer_bins[itv.data][itv.length()] += 1
                print(chrom, itv.begin, itv.end, itv.data, file=args.outfile, sep="\t")
    else:
        df_homopolymers = pl.read_csv(
            args.homopolymers,
            has_header=False,
            separator="\t",
            new_columns=["chrom", "st", "end", "name"],
        )
        homopolymers = defaultdict(it.IntervalTree)
        for chrom, st, end, nt in df_homopolymers.iter_rows():
            all_homopolymer_bins[nt][end - st] += 1
            homopolymers[chrom].add(it.Interval(st, end, nt))

        print("Homopolymers loaded from homopolymers bed.", file=sys.stderr)

    # Multi-histogram colored by nucleotide
    homopolymer_bins = defaultdict(Counter)

    fh_fai = pyfaidx.Faidx(args.reference) if args.reference else None

    for line in df_nucflag_homopolymer_calls.iter_rows(named=True):
        chrom, st, end = line["#chrom"], line["chromStart"], line["chromEnd"]
        homopolymers_chrom = homopolymers.get(chrom)
        if not homopolymers_chrom:
            continue

        ovl = homopolymers_chrom.overlap(st, end)

        assert len(ovl) <= 1, (
            f"Should not have more than one overlap for {chrom}:{st}-{end}."
        )
        # False positive from nucflag (ex. AAAATAAAATAAAAT) or false negative from longdust
        if not ovl:
            len_homopolymer = 0
            if fh_fai:
                rec = fh_fai.fetch(chrom, st + 1, end)
                nt_homopolymer = rec.seq[0]
            else:
                nt_homopolymer = "N"
            if args.outfile_incorrect:
                print(chrom, st, end, sep="\t", file=args.outfile_incorrect)
        else:
            itv_homopolymer: it.Interval = next(iter(ovl))
            len_homopolymer = itv_homopolymer.length()
            nt_homopolymer = itv_homopolymer.data

        homopolymer_bins[nt_homopolymer][len_homopolymer] += 1

    print("Homopolymers binned.", file=sys.stderr)
    # NucFlag overlap
    plot_homopolymers_stacked(
        homopolymer_bins, os.path.join(f"{args.plot_prefix}_ovl_stacked.png")
    )
    plot_homopolymers_split(
        homopolymer_bins,
        os.path.join(f"{args.plot_prefix}_ovl_split.png"),
        include_n=args.reference is None,
    )
    # All
    plot_homopolymers_stacked(
        all_homopolymer_bins,
        os.path.join(f"{args.plot_prefix}_all_stacked.png"),
        ylabel="Number of longdust homopolymer calls",
        # Draw other nucflag overlap as black
        other_homopolymer_bins=homopolymer_bins,
        other_homopolymer_bin_color="black",
        other_homopolymer_bin_label="NucFlag",
    )
    plot_homopolymers_split(
        all_homopolymer_bins,
        os.path.join(f"{args.plot_prefix}_all_split.png"),
        ylabel="Number of longdust homopolymer calls",
        include_n=args.reference is None,
    )

    print(f"Figures saved to {args.plot_prefix}_*.png.", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
