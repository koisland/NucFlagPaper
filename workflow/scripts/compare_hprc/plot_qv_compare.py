import argparse

import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes

COLS_QV = ("chrom", "start", "end", "qv", "bp_err", "bp_correct")
COLS_CTG_MAP = ("chrom_y", "chrom_x", "bp_match")
MARKER_SCALE_FCT = 100_000
MARKER_SIZE_EXAMPLES = [1e6, 5e6, 1e7, 2.5e7]
# https://humanpangenome.org/samples/
AFR_POP_LABELS = {"ACB", "GWD", "ESN", "MKK", "LWK", "ASL", "YRI", "MSL", "ASW"}
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    alignment="left",
)

plt.rcParams["font.family"] = "Arial"


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
    ap.add_argument(
        "-la", "--label_a", type=str, default="HPRC Release 1", help="Label for a."
    )
    ap.add_argument(
        "-lb", "--label_b", type=str, default="HRPC Release 2", help="Label for b."
    )
    ap.add_argument("-ca", "--color_a", type=str, default="red", help="Color for a.")
    ap.add_argument("-cb", "--color_b", type=str, default="blue", help="Color for b.")
    ap.add_argument("-m", "--metadata", type=str, default=None)
    ap.add_argument("-o", "--output_prefix", type=str, default="out")

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
            lbl=pl.lit("a"),
            mtch=pl.col("chrom").str.extract_groups(
                r"^(?<sample>.*?)#(?<hap>.*?)#.*?$"
            ),
        )
        .unnest("mtch")
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
            lbl=pl.lit("b"),
            mtch=pl.col("chrom").str.extract_groups(
                r"^(?<sample>.*?)#(?<hap>.*?)#.*?$"
            ),
        )
        .unnest("mtch")
    )
    labels = {
        "a": args.label_a,
        "b": args.label_b,
    }
    colors = {
        args.label_a: args.color_a,
        args.label_b: args.color_b,
    }

    df_metadata = pl.read_csv(args.metadata, separator="\t", has_header=True)
    df_qvs = pl.concat([df_qv_a, df_qv_b]).join(
        df_metadata, left_on="sample", right_on="Sample ID", how="left"
    )
    df_qvs = df_qvs.with_columns(
        # Set infinite values to median
        pl.when(pl.col("qv").is_infinite())
        .then(pl.col("qv").median().over(["sample", "hap", "lbl"]))
        .otherwise(pl.col("qv"))
        .alias("qv"),
        pl.when(pl.col("Population Abbreviation").is_in(AFR_POP_LABELS))
        .then(pl.lit("AFR"))
        .otherwise(pl.lit("Non-AFR"))
        .alias("Sample"),
        ((pl.col("end") - pl.col("start")) / 1_000_000).alias("Length (Mbp)"),
        # make a numerical column and add some jitter
        # https://stackoverflow.com/a/75541978
        x=pl.when(pl.col("lbl").eq("a")).then(pl.lit(0)).otherwise(pl.lit(1))
        + np.random.uniform(-0.2, 0.2, len(df_qvs)),
        lbl=pl.col("lbl").replace(labels),
    )

    fig, ax = plt.subplots(figsize=(12, 16))
    ax: Axes
    sns.violinplot(
        df_qvs,
        x="lbl",
        y="qv",
        hue="lbl",
        inner="quart",
        palette=colors,
        legend=None,
        ax=ax,
    )
    # https://stackoverflow.com/a/70715319
    # Draw labels on quartiles
    for line in ax.lines:
        data = line.get_data()
        ax.text(
            data[0][data[0].nonzero()][0],
            data[1][0],
            f"{data[1][0]:.0f}",
            size=18,
            zorder=2,
            path_effects=[pe.withStroke(linewidth=3, foreground="w")],
        )

    sns.scatterplot(
        data=df_qvs,
        x="x",
        y="qv",
        hue="Sample",
        size="Length (Mbp)",
        sizes=(1, 800),
        linewidth=1,
        edgecolor="black",
        palette={"AFR": "orange", "Non-AFR": "gray"},
        ax=ax,
    )
    ax.set_xlabel(None)
    ax.set_ylabel("QV (NucFlag)", fontsize=18)
    ax.tick_params(axis="both", which="major", labelsize=18)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    sns.move_legend(
        ax, loc="center", fontsize=18, bbox_to_anchor=(1.125, 0.5), **LEGEND_KWARGS
    )
    df_qvs_summary = df_qvs.group_by(["lbl"]).agg(
        perc_25=pl.col("qv").quantile(0.25),
        median=pl.col("qv").median(),
        mean=pl.col("qv").mean(),
        stdev=pl.col("qv").std(),
        perc_75=pl.col("qv").quantile(0.75),
        mode=pl.col("qv").mode().cast(pl.String).str.join(","),
    )

    df_qvs.write_csv(f"{args.output_prefix}.tsv", separator="\t", include_header=True)
    df_qvs_summary.write_csv(
        f"{args.output_prefix}_summary.tsv", separator="\t", include_header=True
    )
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
