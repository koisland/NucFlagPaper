import sys
import argparse

import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes
from scipy.signal import find_peaks
from scipy.stats import ks_2samp

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


def draw_signif_brackets(
    ax: Axes, x1: int, x2: int, level: float, ylim: tuple[float, float], p: float
):
    # What level is this bar among the bars above the plot?
    # Plot the bar
    y_range = ylim[1] - ylim[0]
    bar_height = (y_range * 0.09 * level) + ylim[1]
    bar_tips = bar_height - (y_range * 0.02)
    ax.plot(
        [x1, x1, x2, x2],
        [bar_tips, bar_height, bar_height, bar_tips],
        lw=1,
        c="k",
        clip_on=False,
    )
    # Significance level
    if p < 0.001:
        sig_symbol = "***"
    elif p < 0.01:
        sig_symbol = "**"
    elif p < 0.05:
        sig_symbol = "*"
    text_height = bar_height + (y_range * 0.01)
    ax.text(
        (x1 + x2) * 0.5,
        text_height,
        sig_symbol,
        ha="center",
        va="bottom",
        c="k",
        fontsize=14,
    )


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
        "-la", "--label_a", type=str, default="Release 1", help="Label for a."
    )
    ap.add_argument(
        "-lb", "--label_b", type=str, default="Release 2", help="Label for b."
    )
    ap.add_argument("-ca", "--color_a", type=str, default="red", help="Color for a.")
    ap.add_argument("-cb", "--color_b", type=str, default="blue", help="Color for b.")
    ap.add_argument("-m", "--metadata", type=str, default=None)
    ap.add_argument("-o", "--output_prefix", type=str, default="out")
    ap.add_argument(
        "-t", "--type", choices=["length", "qv"], type=str, default="length"
    )

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
        0: args.label_a,
        1: args.label_b,
    }
    colors = {
        args.label_a: args.color_a,
        args.label_b: args.color_b,
    }
    df_metadata = pl.read_csv(args.metadata, separator="\t", has_header=True)
    df_qvs = pl.concat([df_qv_a, df_qv_b]).join(
        df_metadata, left_on="sample", right_on="Sample ID", how="left"
    )
    if args.type == "length":
        col = "Length (Mbp)"
        ylabel = "Contig length (Mbp)"
    else:
        col = "qv"
        ylabel = "QV (NucFlag)"
    df_qvs = (
        df_qvs.with_columns(
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
        )
        .with_columns(
            group=pl.when(pl.col("lbl").eq("a")).then(pl.lit(0)).otherwise(pl.lit(1))
        )
        .with_columns(
            Group=pl.col("group").cast(pl.String).replace(labels),
        )
    )

    g = sns.catplot(
        data=df_qvs,
        x="Group",
        y=col,
        hue="Group",
        hue_order=colors.keys(),
        order=colors.keys(),
        inner="quart",
        palette=colors,
        legend=None,
        kind="violin",
        height=8,
        aspect=0.5,
        alpha=0.7,
    )
    ax: Axes = g.ax
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

    # One-sided test two-sample KS test
    r1_vals, r2_vals = df_qvs.select("Group", col).partition_by("Group")
    res = ks_2samp(r1_vals[col], r2_vals[col], alternative="greater")
    draw_signif_brackets(ax, x1=0, x2=1, level=1.0, ylim=ax.get_ylim(), p=res.pvalue)
    print(f"R2 has greater {args.type} than R1 ({res.pvalue})", file=sys.stderr)
    sns.stripplot(
        data=df_qvs,
        x="Group",
        y=col,
        hue="Group",
        palette=colors,
        jitter=False,
        linewidth=0.5,
        size=4,
        edgecolor="black",
        ax=ax,
        zorder=1,
    )
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.tick_params(axis="both", which="major", labelsize=18)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    df_qvs = df_qvs.with_columns(Group=pl.col("Group").str.replace("\n", " "))

    peak_rows = []
    for grp, df_grp in df_qvs.group_by(["Group"]):
        counts = df_grp[col].round().value_counts().sort(by=col)
        peaks, props = find_peaks(counts["count"], height=10)
        for pk, ht in zip(peaks, props["peak_heights"]):
            peak_rows.append((pk, ht, grp[0]))

    df_peaks = pl.DataFrame(
        peak_rows,
        orient="row",
        schema={col: pl.Float32, "ht": pl.Int32, "Group": pl.String},
        infer_schema_length=None,
    )
    df_peaks.write_csv(
        f"{args.output_prefix}_peaks.tsv", separator="\t", include_header=True
    )

    # sns.move_legend(
    #     ax, loc="center", fontsize=18, bbox_to_anchor=(1.125, 0.5), title=None, **LEGEND_KWARGS
    # )
    df_qvs_summary = df_qvs.group_by(["Group"]).agg(
        perc_25=pl.col(col).quantile(0.25),
        median=pl.col(col).median(),
        mean=pl.col(col).mean(),
        stdev=pl.col(col).std(),
        perc_75=pl.col(col).quantile(0.75),
    )

    df_qvs.write_csv(f"{args.output_prefix}.tsv", separator="\t", include_header=True)
    df_qvs_summary.write_csv(
        f"{args.output_prefix}_summary.tsv", separator="\t", include_header=True
    )
    g.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
