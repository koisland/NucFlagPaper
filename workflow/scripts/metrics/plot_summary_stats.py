import sys
import ast
import argparse

import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes

OUTFILE_KWARGS = dict(file=sys.stdout, sep="\t")
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    title=None,
    frameon=False,
)


# https://www.statology.org/what-is-a-good-f1-score/
# F1 Score = 2 * (Precision * Recall) / (Precision + Recall)
def calculate_metrics(df: pl.DataFrame) -> tuple[float, float, float]:
    prec = df["precision"].mean()
    recall = df["recall"].mean()
    f1 = 2 * (prec * recall) / (prec + recall)
    return prec, recall, f1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i", "--input_summary_globs", help="Input summary globs", nargs="+"
    )
    ap.add_argument("-l", "--labels", help="Labels", nargs="+")
    ap.add_argument("-c", "--colors", help="Colors", nargs="+")
    ap.add_argument("-s", "--figsize", help="Figure size.", type=str, default="(16, 8)")
    ap.add_argument("-o", "--output", help="Output plot", default="out.png")
    args = ap.parse_args()

    labels = []
    metrics = []
    types = []
    print("label", "precision", "recall", "f1", **OUTFILE_KWARGS)
    for glob_file, label in zip(args.input_summary_globs, args.labels, strict=True):
        df_summaries = pl.read_csv(
            glob_file,
            separator="\t",
            glob=True,
            has_header=True,
        )
        prec, recall, f1 = calculate_metrics(df_summaries)

        labels.extend((label, label, label))
        metrics.extend((prec, recall, f1))
        types.extend(("Precision", "Recall", "F1"))

        print(label, prec, recall, f1, **OUTFILE_KWARGS)

    df = pl.DataFrame({"label": labels, "percent": metrics, "type": types})
    fig, axes = plt.subplots(
        nrows=1, ncols=3, layout="tight", figsize=ast.literal_eval(args.figsize)
    )
    label_colors = dict(zip(args.labels, args.colors, strict=True))

    for i, col in enumerate(("Precision", "Recall", "F1")):
        ax: Axes = axes[i]
        df_type = df.filter(pl.col("type") == col)

        sns.barplot(
            df_type,
            x="label",
            y="percent",
            hue="label",
            palette=label_colors,
            legend=None,
            ax=ax,
        )
        for bar in ax.containers:
            ax.bar_label(bar, fontsize=8, fmt=lambda x: f"{x * 100:.1f}%")

        yticks = [tick / 100.0 for tick in range(0, 120, 20)]
        ax.set_yticks(yticks, [f"{tick * 100.0:.1f}" for tick in yticks])

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        for lbl in ax.xaxis.get_majorticklabels():
            lbl.set_rotation(45)
            lbl.set_horizontalalignment("right")
            lbl.set_rotation_mode("anchor")
            lbl.set_color(label_colors.get(lbl.get_text(), "white"))
            lbl.set_path_effects([pe.Stroke(linewidth=0.2, foreground="black")])

        ax.set_xlabel(None)
        ax.set_ylabel("Percent (%)")
        ax.set_title(col, fontsize="xx-large")

    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
