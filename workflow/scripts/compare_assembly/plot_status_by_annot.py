import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes
from matplotlib.patches import Patch

plt.rcParams["font.family"] = "Arial"


LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)


def format_bar_ax(ax: Axes, name_colors: dict[str, str]):
    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)

    for lbl in ax.xaxis.get_majorticklabels():
        lbl_text = lbl.get_text()
        lbl.set_path_effects([pe.Stroke(linewidth=0.1, foreground="black")])
        lbl.set_color(name_colors[lbl_text])
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")

    for c in ax.containers:
        labels = [str(round(v.get_height())) for v in c]
        ax.bar_label(c, labels=labels, label_type="edge", fontsize=9)

    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_ylabel("QV (NucFlag)", fontsize=14)
    ax.set_xlabel(None)


def format_group_length_ax(ax: Axes, name_colors: dict[str, str]):
    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)

    ax.set_ylabel("Length (Mbp)", fontsize=14)
    ax.set_xlabel(None)
    ax.tick_params(axis="both", which="major", labelsize=14)

    for lbl in ax.xaxis.get_majorticklabels():
        lbl_text = lbl.get_text()
        lbl.set_path_effects([pe.Stroke(linewidth=0.1, foreground="black")])
        lbl.set_color(name_colors[lbl_text])
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")

    for c in ax.containers:
        labels = [str(round(v.get_height())) for v in c]
        ax.bar_label(c, labels=labels, label_type="edge", fontsize=9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i", "--input_statuses", type=str, nargs="+", help="Input NucFlag statuses"
    )
    ap.add_argument("-l", "--labels", type=str, nargs="+", help="Labels.")
    ap.add_argument("-c", "--colors", type=str, nargs="+", help="Colors for labels")
    ap.add_argument(
        "-g",
        "--groups",
        type=str,
        default=None,
        help="TSV with groups. Order by position and color in 2nd column.",
    )
    ap.add_argument("-t", "--title", type=str, default=None, help="Output plot title")
    ap.add_argument(
        "-o", "--output_prefix", type=str, default="out", help="Output plot prefix"
    )
    args = ap.parse_args()

    dfs = []
    for label, status in zip(args.labels, args.input_statuses, strict=True):
        df = pl.read_csv(status, has_header=True, separator="\t").with_columns(
            label=pl.lit(label)
        )
        dfs.append(df)

    label_colors = dict(zip(args.labels, args.colors, strict=True))
    df_all = pl.concat(dfs)
    if args.groups:
        df_groups = pl.read_csv(
            args.groups,
            separator="\t",
            has_header=False,
            new_columns=["group", "color"],
        )
        group_colors: dict[str, str] = dict(df_groups.iter_rows())
    else:
        all_groups = df_all["group"].unique().sort().to_list()
        group_colors = dict(zip(all_groups, ["black"] * len(all_groups)))

    width_scale_factor = 1
    width = width_scale_factor * len(group_colors)
    fig, axes = plt.subplots(
        figsize=(width, 6), layout="constrained", nrows=2, ncols=1, sharex=True
    )

    # Get intersection between
    all_groups = set(df_all["group"].unique())
    group_order = [grp for grp in group_colors.keys() if grp in all_groups]

    ax_bar: Axes = axes[0]
    sns.barplot(
        df_all,
        x="group",
        hue="label",
        y="QV",
        ax=ax_bar,
        legend=None,
        order=group_order,
        hue_order=label_colors.keys(),
        palette=label_colors,
    )
    format_bar_ax(ax_bar, group_colors)

    ax_liftover: Axes = axes[1]
    sns.barplot(
        df_all.with_columns(pl.col("group_length") / 1_000_000),
        x="group",
        hue="label",
        y="group_length",
        palette=label_colors,
        order=group_order,
        hue_order=label_colors.keys(),
        ax=ax_liftover,
        legend=None,
    )
    format_group_length_ax(ax_liftover, name_colors=group_colors)

    fig_legend = plt.figure(layout="constrained", figsize=(2, 2))
    fig_legend.legend(
        handles=[
            Patch(edgecolor="black", facecolor=color) for color in label_colors.values()
        ],
        labels=label_colors.keys(),
        fontsize=14,
        loc="center",
        **LEGEND_KWARGS,
    )
    fig.suptitle(
        args.title,
        x=0.0,
        fontsize=18,
        weight="bold",
        horizontalalignment="left",
        verticalalignment="center",
    )
    df_all.write_csv(f"{args.output_prefix}.tsv", separator="\t")
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=300)
    fig_legend.savefig(f"{args.output_prefix}_legend.png", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
