import bisect
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from collections import Counter

MAPQ_COLORS = {
    "0": "#ffffff",
    "1-5": "#666666",
    "5-10": "#8f59a7",
    "10-15": "#5954a8",
    "15-20": "#01aef3",
    "20-25": "#04b99e",
    "25-30": "#8bc83b",
    "30-35": "#cdde25",
    "35-40": "#fff600",
    "40-45": "#ffc309",
    "45-50": "#fa931a",
    "50-55": "#f8631f",
    "55-60": "#f21821",
    "60": "#FF8DA1",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input calls and mapq bedgraph intersected (-wa -wb). BED9+BED4",
    )
    ap.add_argument("-o", "--output", required=True, help="Outfile.")
    args = ap.parse_args()

    df_calls = pl.read_csv(
        args.input,
        separator="\t",
        has_header=False,
        new_columns=(
            "chrom",
            "st",
            "end",
            "name",
            "strand",
            "score",
            "tst",
            "tend",
            "item_rgb",
            "lchrom",
            "lst",
            "lend",
            "mapq",
        ),
    ).filter(~pl.col("chrom").is_in(["chrEBV", "chrM"]))

    mapq_counter = Counter()
    breakpoints = list(range(0, 70, 5))
    for grp, df_grp in df_calls.group_by(["chrom", "st", "end", "name"]):
        chrom, st, end, name = grp
        mapq = df_grp["mapq"].min()

        if mapq == 60.0:
            mapq_label = "60"
        elif mapq == 0.0:
            mapq_label = "0"
        else:
            idx = bisect.bisect(breakpoints, mapq)
            mapq_lower_bound = breakpoints[idx - 1]
            # Unreachable due to above. Set to 1 so descriptive
            if mapq_lower_bound == 0:
                mapq_lower_bound = 1
            mapq_upper_bound = breakpoints[idx]
            mapq_label = f"{mapq_lower_bound}-{mapq_upper_bound}"

        mapq_counter[mapq_label] += 1

    df_mapq_counts = (
        pl.DataFrame(
            [
                {"label": lbl, "count": cnt, "first": int(lbl.split("-")[0])}
                for lbl, cnt in mapq_counter.items()
            ],
            orient="row",
        )
        .sort(by="first")
        .drop("first")
    )
    fig, ax = plt.subplots(figsize=(16, 8), layout="constrained")
    sns.barplot(
        x=df_mapq_counts["label"],
        y=df_mapq_counts["count"],
        hue=df_mapq_counts["label"],
        palette=MAPQ_COLORS,
        linewidth=1.0,
        edgecolor="black",
        legend=False,
        ax=ax,
    )

    for cnter in ax.containers:
        ax.bar_label(
            cnter,
            label_type="center",
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        )

    ax.set_xlabel("MAPQ range for Element BAM")
    ax.set_ylabel("Number of intersected calls")
    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)

    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
