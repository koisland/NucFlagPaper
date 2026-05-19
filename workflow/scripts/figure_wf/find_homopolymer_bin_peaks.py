import os
import sys
import polars as pl
import scipy.signal as sig
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def draw_fig_with_peaks(
    ax: Axes, df: pl.DataFrame, height: int = 1000, prom: int = 5
) -> None:
    peaks = sig.find_peaks(df["count"], prominence=prom, height=height)

    ax.bar(df["len"], df["count"])
    for peak, ht in zip(df["len"][peaks[0]], peaks[1]["peak_heights"]):
        ax.axvline(peak, linestyle="dotted", color="black")
        ax.annotate(f"{peak} ({ht})", xy=(peak, ht))

    if len(peaks[0]) == 0:
        length, count = df.filter(pl.col("count").eq(pl.col("count").max())).row()
        ax.axvline(length, linestyle="dotted", color="black")
        ax.annotate(f"{length} ({count})", xy=(length, count))

    df_non_zero = df.filter(pl.col("len") != 0)
    prop_gt_20 = (
        df_non_zero.filter(pl.col("len") >= 20)["count"].sum()
        / df_non_zero["count"].sum()
    )
    ax.set_title(f"Proportion >= 20 bp ({prop_gt_20})")


NT_ORDER = ("A", "T", "G", "C")
HEIGHTS = {"A": 500, "T": 500, "G": 20, "C": 20}
PROMS = {"A": 500, "T": 500, "G": 20, "C": 20}


def main():
    # Homopolymer bins TSV (nt, length, count)
    bins = sys.argv[1]
    wd = os.path.dirname(bins)
    fname = os.path.splitext(os.path.basename(bins))[0]
    df = pl.read_csv(bins, separator="\t", new_columns=["nt", "len", "count"])
    df_agg = df.group_by(["len"]).agg(pl.col("count").sum()).sort("len")

    fig_all, ax_all = plt.subplots(figsize=(6, 4), layout="tight")
    draw_fig_with_peaks(ax_all, df_agg)

    fig_split, axes_split = plt.subplots(
        ncols=4, nrows=1, figsize=(16, 4), layout="tight", sharey=True, sharex=True
    )
    axes_split: list[Axes]

    for i, ax in enumerate(axes_split):
        nt = NT_ORDER[i]
        df_nt = df.filter(pl.col("nt").eq(pl.lit(nt)))
        df_nt_agg = df_nt.group_by(["len"]).agg(pl.col("count").sum()).sort("len")
        draw_fig_with_peaks(ax, df_nt_agg, height=HEIGHTS[nt], prom=PROMS[nt])
        ax.set_title(f"{nt}\n{ax.get_title()}")

    fig_all.savefig(os.path.join(wd, f"{fname}.png"), bbox_inches="tight")
    fig_split.savefig(os.path.join(wd, f"{fname}_split.png"), bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
