import argparse
import sys
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from collections import deque
from itertools import islice
from bisect import bisect


# https://docs.python.org/3/library/itertools.html
def sliding_window(iterable, n):
    "Collect data into overlapping fixed-length chunks or blocks."
    # sliding_window('ABCDEFG', 3) → ABC BCD CDE DEF EFG
    iterator = iter(iterable)
    window = deque(islice(iterator, n - 1), maxlen=n)
    for x in iterator:
        window.append(x)
        yield tuple(window)


# Format from https://www.researchgate.net/figure/Assembly-results-for-four-assemblers-and-three-human-samples-before-polishing-a-NGx_fig2_334713255
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_fai", type=str, nargs="+")
    ap.add_argument("-l", "--labels", type=str, nargs="+")
    ap.add_argument("-c", "--colors", type=str, nargs="+")
    ap.add_argument("-o", "--output", type=str, default="out.png")
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
    max_y = 0
    for label, file, color in zip(
        args.labels, args.input_fai, args.colors, strict=True
    ):
        df = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
                columns=[0, 1],
                new_columns=["contig", "contig_length"],
            )
            .sort(by="contig_length", descending=True)
            .with_columns(
                contig_length=pl.col("contig_length") / 1_000_000,
                coverage=(
                    pl.col("contig_length").cum_sum() / pl.col("contig_length").sum()
                )
                * 100,
            )
        )

        if color == "None":
            color = "black"

        # x - cumulative coverage, y - contig length
        x = 0
        for rows in sliding_window(iterable=df.iter_rows(), n=2):
            row_1, row_2 = rows
            _, length_1, cov_1 = row_1
            _, length_2, _ = row_2
            # ---+
            x_pts_1 = [x, cov_1]
            y_pts_1 = [length_1, length_1]
            ax.plot(x_pts_1, y_pts_1, label=label, color=color)
            # ---+
            #    |
            x_pts_2 = [cov_1, cov_1]
            y_pts_2 = [length_1, length_2]
            ax.plot(x_pts_2, y_pts_2, label=label, color=color)
            x = cov_1
            max_y = max(length_1, max_y)

        # https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
        aun = (df["contig_length"] ** 2).sum() / df["contig_length"].sum()
        print(f"{label} auN: {aun}", file=sys.stderr)
        # Draw ng50 label
        idx = bisect(df["coverage"], 50.0)
        _, y_ng50, x_ng50_cov = df.row(idx)
        ax.annotate(
            f"NG50: {y_ng50:.1f} Mbp",
            (x_ng50_cov, y_ng50),
            color=color,
            path_effects=[pe.withStroke(linewidth=4.0, foreground="white")],
        )
        # ax.axhline(y=y_ng50, linestyle="dotted", color=color)
        # ax.axhline(y=aun, linestyle="dashed", color=color)

    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(0.0, max_y + 5.0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    labels_handles = dict(zip(labels, handles))
    ax.legend(
        labels=labels_handles.keys(),
        handles=labels_handles.values(),
        loc="upper right",
        fancybox=False,
        frameon=False,
    )
    ax.set_xlabel(r"% of genome")
    ax.set_ylabel("Contig length (Mbp)")
    fig.savefig(args.output, bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
