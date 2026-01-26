import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.patches import Patch
from matplotlib.colors import rgb2hex

COLS_CALLS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "tstart",
    "tend",
    "item_rgb",
)
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    title=None,
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    fancybox=False,
    frameon=False,
)
DIPLOID_GENOME_SIZE = 6_200_000_000


def draw_combined_err_diff(df: pl.DataFrame, outfile: str, name_colors: dict[str, str]):
    name_fig, ax_fig = plt.subplots(layout="constrained", figsize=(8, 8))
    name_bars = ax_fig.bar(
        x=df["name"],
        height=df["diff"],
        facecolor=[name_colors[name] for name in df["name"]],
        edgecolor="black",
    )
    ax_fig.bar_label(
        name_bars,
        fontsize=8,
        fmt=lambda length: f"{length / 1_000_000:.1f}",
        label_type="center",
        path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        padding=3,
    )
    ax_fig.yaxis.set_major_formatter(lambda x, _: f"{x / 1_000_000:.1f}")
    ax_fig.tick_params(axis="x", labelrotation=45)
    ax_fig.margins(x=0.01)
    ax_fig.set_ylabel("Combined length (Mbp)")
    for spine in ("top", "right"):
        ax_fig.spines[spine].set_visible(False)

    for label in ax_fig.get_xticklabels():
        label.set_color(name_colors[label.get_text()])
        label.set_path_effects([pe.withStroke(linewidth=0.2, foreground="black")])

    name_fig.savefig(outfile, bbox_inches="tight", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-a",
        "--input_a",
        type=str,
        help="Input calls for group a. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument(
        "-b",
        "--input_b",
        type=str,
        help="Input calls for group b. Expects PanSN naming spec with # as delimiter.",
    )
    ap.add_argument("-o", "--output", default="output.png")
    ap.add_argument("-n", "--output_name", default="output_name.png")
    ap.add_argument("-nd", "--output_name_diff", default="out_name_diff.tsv")
    ap.add_argument("-sd", "--output_sm_diff", default="out_sm_diff.tsv")
    ap.add_argument("-snd", "--output_sm_name_diff", default="out_sm_name_diff.tsv")

    args = ap.parse_args()

    df_a = (
        pl.read_csv(
            args.input_a,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_CALLS,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("a"), length=pl.col("end") - pl.col("start"))
        .filter(pl.col("name").ne(pl.lit("correct")))
    )
    df_b = (
        pl.read_csv(
            args.input_b,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=COLS_CALLS,
        )
        .with_columns(
            mtch=pl.col("chrom").str.extract_groups(r"^(?<sample>.*?)#(?<hap>.*?)#.*?$")
        )
        .unnest("mtch")
        .with_columns(label=pl.lit("b"), length=pl.col("end") - pl.col("start"))
        .filter(pl.col("name").ne(pl.lit("correct")))
    )

    df_all = pl.concat((df_a, df_b))
    df_grp_no_hap = df_all.group_by(["sample", "label", "name"]).agg(
        pl.col("length").sum()
    )
    df_grp_no_hap_wide = (
        df_grp_no_hap.pivot(on="name", index=("sample", "label"), values="length")
        .sort(by=("sample", "label"))
        .fill_null(0)
    )

    unique_samples = df_grp_no_hap_wide["sample"].unique().sort()
    name_colors = {
        name: rgb2hex(tuple(int(c) / 255 for c in color.split(",")))
        for name, color in df_all.select("name", "item_rgb").iter_rows()
    }
    label_samples = np.arange(len(unique_samples))
    bottoms = {
        "a": np.zeros(len(unique_samples)),
        "b": np.zeros(len(unique_samples)),
    }

    fig, ax = plt.subplots(layout="constrained", figsize=(24, 6))
    width = 0.3
    multiplier_adj = 1.1

    for call_name, call_color in name_colors.items():
        multiplier = 0
        for label in bottoms.keys():
            bottom = bottoms[label]
            offset = width * multiplier
            lengths = df_grp_no_hap_wide.filter(pl.col("label").eq(label))[
                call_name
            ].to_numpy()
            rects = ax.bar(
                label_samples + offset,
                lengths,
                width,
                label=call_name,
                bottom=bottom,
                facecolor=call_color,
                edgecolor="black",
            )
            # Don't plot small regions which you wouldn't see anyways.
            if lengths.mean() > 10_000:
                ax.bar_label(
                    rects,
                    fmt=lambda length: f"{length / 1_000_000:.1f}",
                    fontsize=8,
                    label_type="center",
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
                    padding=3,
                )

            # Label colors should be red/blue for r1/r2 based on palette
            multiplier += multiplier_adj
            bottom += lengths

    # Minimize spacing between xtick and spine
    ax.margins(x=0.01)

    # Backwards since last one is first one in
    ax.legend(
        labels=reversed(name_colors.keys()),
        handles=[
            Patch(edgecolor="black", facecolor=color)
            for color in reversed(name_colors.values())
        ],
        **LEGEND_KWARGS,
    )
    # Nudge xtick slightly so centered.
    ax.set_xticks(label_samples + 0.15, unique_samples)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.tick_params(axis="x", labelrotation=45)

    # Just eyeballing
    ax.yaxis.set_major_formatter(lambda x, _: f"{x / 1_000_000:.1f}")
    ax.set_ylabel("Cumulative length (Mbp)")
    fig.savefig(args.output, bbox_inches="tight", dpi=150)

    # Agg difference in called bases by sample and sample/name
    df_cmp_ab = (
        df_grp_no_hap.group_by(["sample", "name", "label"])
        .agg(pl.col("length").sum())
        .pivot(on="label", index=["sample", "name"])
        .select("sample", "name", "a", "b")
        .with_columns(diff=pl.col("b") - pl.col("a"))
    )
    df_cmp_name_ab = (
        df_cmp_ab.group_by(["name"])
        .agg(a=pl.col("a").sum(), b=pl.col("b").sum())
        .with_columns(diff=pl.col("b") - pl.col("a"))
        .sort(by="name")
    )
    df_cmp_sm_ab = (
        df_cmp_ab.group_by(["sample"])
        .agg(a=pl.col("a").sum(), b=pl.col("b").sum())
        .with_columns(diff=pl.col("b") - pl.col("a"))
    )

    draw_combined_err_diff(
        df_cmp_name_ab, outfile=args.output_name, name_colors=name_colors
    )
    median_diff_called_bases = df_cmp_sm_ab["diff"].median()
    perc_median_diff_called_bases = (
        median_diff_called_bases / DIPLOID_GENOME_SIZE
    ) * 100
    print(
        f"Median reduction in called bases (Mbp) from a to b: {median_diff_called_bases} ({perc_median_diff_called_bases}%)"
    )
    df_cmp_ab.write_csv(args.output_sm_name_diff, separator="\t")
    df_cmp_sm_ab.write_csv(args.output_sm_diff, separator="\t")
    df_cmp_name_ab.write_csv(args.output_name_diff, separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
