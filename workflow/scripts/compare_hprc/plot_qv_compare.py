import argparse

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D


COLS_QV = ("chrom", "start", "end", "qv", "bp_err", "bp_correct")
COLS_CTG_MAP = ("chrom_y", "chrom_x", "bp_match")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-a",
        "--qv_a",
        type=str,
        help="QV for group a. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument(
        "-b",
        "--qv_b",
        type=str,
        help="QV for group b. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument("-ca", "--color_a", type=str, default="red", help="Color for a.")
    ap.add_argument("-cb", "--color_b", type=str, default="blue", help="Color for b.")
    ap.add_argument("-o", "--output", type=str, default="out.png")

    args = ap.parse_args()
    df_qv_a = (
        pl.read_csv(
            args.qv_a,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_QV,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("a"))
    )
    df_qv_b = (
        pl.read_csv(
            args.qv_b,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_QV,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("b"))
    )
    df_qv_all = pl.concat([df_qv_a, df_qv_b])
    colors = {
        "a": args.color_a,
        "b": args.color_b,
    }
    fig, ax = plt.subplots(figsize=(8, 8))
    for sm, df_grp in df_qv_all.group_by(["sample"]):
        df_grp_agg = (
            df_grp.with_columns(
                pl.when(pl.col("qv").is_infinite())
                .then(pl.col("qv").median().over(["hap", "label"]))
                .otherwise(pl.col("qv"))
            )
            .group_by(["hap", "label"])
            .agg(
                qv=pl.col("qv").mean(),
                len=(pl.col("end") - pl.col("start")).mean().cast(pl.UInt64),
            )
            .sort("label", "hap")
        )
        df_points = (
            df_grp_agg.pivot(on="hap", index="label", values=["qv", "len"])
            .with_columns(len=(pl.col("len_1") + pl.col("len_2")) / 2)
            .rename({"qv_1": "x", "qv_2": "y"})
            .select("label", "x", "y", "len")
        )
        for point in df_points.iter_rows(named=True):
            size = point["len"] // 100_000
            color = colors[point["label"]]
            ax.scatter(
                x=point["x"],
                y=point["y"],
                s=size,
                facecolor="none",
                edgecolor=color,
                linestyle="--",
                marker="o",
            )

        # Arrow from a -> b.
        row_a = df_points.filter(pl.col("label") == "a").row(0, named=True)
        row_b = df_points.filter(pl.col("label") == "b").row(0, named=True)
        arr = FancyArrowPatch(
            (row_a["x"], row_a["y"]),
            (row_b["x"], row_b["y"]),
            arrowstyle="-|>",
            color="black",
            mutation_scale=5,
        )
        ax.add_patch(arr)

    # TODO: Add legend with bubble scaling
    ax.legend(
        labels=[
            "HPRC Release 1",
            "HPRC Release 2",
        ],
        handles=[
            Line2D(
                [],
                [],
                color="white",
                marker="o",
                markerfacecolor="none",
                markeredgecolor="red",
            ),
            Line2D(
                [],
                [],
                color="white",
                marker="o",
                markerfacecolor="none",
                markeredgecolor="blue",
            ),
        ],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
    )
    ax.set_xlabel("QV (Hap1)")
    ax.set_ylabel("QV (Hap2)")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.savefig(args.output, bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
