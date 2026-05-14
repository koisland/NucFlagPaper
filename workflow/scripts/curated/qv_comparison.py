import sys
import argparse
import scipy.stats

from matplotlib.axes import Axes
from matplotlib.lines import Line2D

import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

# https://kyrrego.github.io/blog/2025/01/31/my-first-blog.html
plt.rcParams["font.family"] = "Arial"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-l",
        "--labels",
        nargs="+",
        required=True,
        help="Labels of equal number of x and y.",
        type=str,
    )
    ap.add_argument(
        "-c",
        "--colors",
        nargs="+",
        required=True,
        help="Colors of equal number of x and y.",
        type=str,
    )
    ap.add_argument(
        "-x",
        "--qv_x",
        nargs="+",
        required=True,
        help="QV of x.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-y",
        "--qv_y",
        nargs="+",
        required=True,
        help="QV of y.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument("--xlabel", default="QV (Merqury)", help="x-label", type=str)
    ap.add_argument("--ylabel", default="QV (NucFlag)", help="y-label", type=str)
    ap.add_argument(
        "--highlight-acrocentrics",
        dest="highlight_acrocentrics",
        action="store_true",
        help="Highlight acrocentric chromosomes.",
    )
    ap.add_argument(
        "-o", "--output_prefix", required=True, help="Output plot prefix.", type=str
    )

    # x - nucflag, y - merqury
    args = ap.parse_args()

    def add_regression_line_text(ax: Axes, data: pl.DataFrame, **kws):
        r, p = scipy.stats.pearsonr(data["qv_x"], data["qv_y"])
        ax.text(
            0.05,
            0.8,
            "r={:.2f}, p={:.2g}".format(r, p),
            fontsize="large",
            transform=ax.transAxes,
        )

    for label, color, qv_x, qv_y in zip(
        args.labels, args.colors, args.qv_x, args.qv_y, strict=True
    ):
        # Expects chrom, qv
        df_x = pl.read_csv(
            qv_x,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_x"],
            comment_prefix="#",
            has_header=False,
        )
        df_y = pl.read_csv(
            qv_y,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_y"],
            comment_prefix="#",
            has_header=False,
        )
        # If main chromosomes, plot.
        df_all = df_x.join(df_y, on="chrom", how="inner").filter(
            pl.col("chrom").str.contains_any(["MATERNAL", "PATERNAL"])
        )

        p = sns.jointplot(
            df_all,
            x="qv_x",
            y="qv_y",
            kind="reg",
            xlim=(45, 55),
            ylim=(30, 50),
            marker="o",
            color=color,
            line_kws=dict(color=color, linestyle="dotted"),
        )
        # Highlight how acrocentric chromosomes are the outliers
        if args.highlight_acrocentrics:
            df_acro_chrs = df_all.filter(
                pl.col("chrom").str.contains_any(
                    ["chr13", "chr14", "chr15", "chr21", "chr22"]
                )
            )

            p.ax_joint.scatter(
                df_acro_chrs["qv_x"],
                df_acro_chrs["qv_y"],
                marker="o",
                facecolors="none",
                edgecolors="black",
            )

        # set median.
        median_x = df_x["qv_x"].median()
        median_y = df_y["qv_y"].median()

        ax: Axes = plt.gca()
        ax.axvline(x=median_x, linestyle="dotted", color=color)
        ax.axhline(y=median_y, linestyle="dotted", color=color)
        print(label, median_x, median_y, file=sys.stderr, sep="\t")

        add_regression_line_text(ax, df_all)
        ax.set_title(
            label,
            color=color,
            x=0.05,
            y=0.95,
            ha="left",
            weight="bold",
            fontsize="xx-large",
        )
        ax.set_xlabel(args.xlabel, fontsize="x-large")
        ax.set_ylabel(args.ylabel, fontsize="x-large")
        # annotate(ax, df_other_chrs)

        if args.highlight_acrocentrics:
            ax.legend(
                labels=["Acrocentric"],
                handles=[
                    Line2D(
                        [],
                        [],
                        color="white",
                        marker="o",
                        markerfacecolor="none",
                        markeredgecolor="black",
                    )
                ],
                loc="upper right",
            )
        p.savefig(f"{args.output_prefix}_{label}.pdf", dpi=300)
        p.savefig(f"{args.output_prefix}_{label}.png", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
