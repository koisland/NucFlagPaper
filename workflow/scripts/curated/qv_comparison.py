import scipy.stats
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt


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
    ap.add_argument("--xlabel", default="QV (x)", help="x-label", type=str)
    ap.add_argument("--ylabel", default="QV (y)", help="y-label", type=str)
    ap.add_argument("-o", "--output", required=True, help="Output plot.", type=str)

    # x - nucflag, y - merqury
    args = ap.parse_args()

    fig, ax = plt.subplots(layout="constrained", figsize=(8, 8))
    for i, (label, color, qv_x, qv_y) in enumerate(
        zip(args.labels, args.colors, args.qv_x, args.qv_y, strict=True), 1
    ):
        # Expects chrom, qv
        df_x = pl.read_csv(
            qv_x,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_x"],
            has_header=False,
        )
        df_y = pl.read_csv(
            qv_y,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_y"],
            has_header=False,
        )
        # If inf, set to 100.0
        df_all = df_x.join(df_y, on="chrom", how="inner").with_columns(
            qv_x=pl.when(pl.col("qv_x").is_infinite())
            .then(pl.lit(100.0))
            .otherwise(pl.col("qv_x")),
            qv_y=pl.when(pl.col("qv_y").is_infinite())
            .then(pl.lit(100.0))
            .otherwise(pl.col("qv_y")),
        )

        p = sns.regplot(
            df_all,
            x="qv_x",
            y="qv_y",
            ci=99,
            marker="o",
            color=color,
            label=label,
            ax=ax,
            line_kws=dict(color="black", linestyle="dotted"),
        )
        # https://www.statology.org/seaborn-regplot-equation/
        slope, intercept, r, p, sterr = scipy.stats.linregress(
            x=p.get_lines()[0].get_xdata(), y=p.get_lines()[0].get_ydata()
        )
        # Midpt-ish
        # Offset if multiple
        x = 50 + 11 * i
        y = intercept + (slope * x)
        ax.text(
            x,
            y,
            f"y = {round(intercept, 3)} + {round(slope, 3)}x (r = {r:.3f})",
            color=color,
        )

    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)

    ax.legend(loc="lower right")

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.savefig(args.output, dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
