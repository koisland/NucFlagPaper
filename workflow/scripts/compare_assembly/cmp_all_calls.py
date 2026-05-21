import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex

plt.rcParams["font.family"] = "Arial"

LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    title=None,
)


# https://stackoverflow.com/a/67594395
def update_bars(ax: Axes, round_to: int):
    for c in ax.containers:
        labels = [f"{round(v.get_height() / 1_000_000, round_to)}" for v in c]
        ax.bar_label(c, labels=labels, label_type="edge", fontsize=18)
    ax.set_xlabel(None)


def rgb_to_hex(srs: pl.Series) -> pl.Series:
    color_hex = []
    for elem in srs:
        if elem.startswith("#"):
            color_hex.append(elem)
        else:
            rgb = tuple(int(e) / 255 for e in elem.split(","))
            assert len(rgb) == 3, f"Invalid item_rgb format for {rgb}"
            color_hex.append(rgb2hex(rgb))
    return pl.Series(name="color", values=color_hex)


def draw_bar_all(
    df_all: pl.DataFrame,
    label_colors: dict[str, str],
    name_colors: dict[str, str],
    output_prefix: str,
):
    fig, ax = plt.subplots(
        figsize=(24, 6),
        layout="constrained",
        sharex=True,
        sharey=True,
    )
    ax: Axes
    max_length = df_all["length"].max()
    yticks = list(range(0, int(max_length), 1_000_000))
    yticklabels = [f"{tick / 1_000_000}" for tick in yticks]
    sns.barplot(
        data=df_all,
        x="name",
        y="length",
        hue="lbl",
        order=df_all["name"].sort().to_list(),
        hue_order=label_colors.keys(),
        legend="full",
        palette=label_colors,
        ax=ax,
    )

    update_bars(ax, round_to=1)
    for lbl in ax.xaxis.get_majorticklabels():
        # Color and stroke around text
        lbl.set_color(name_colors[lbl.get_text()])
        lbl.set_fontsize(18)
        lbl.set_path_effects([pe.Stroke(linewidth=0.1, foreground="black")])
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")

    ax.set_ylabel("Length (Mbp)", fontsize=18)
    ax.set_yticks(yticks, yticklabels)
    ax.tick_params(axis="both", which="major", labelsize=18)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    sns.move_legend(ax, loc="upper right", fontsize=18, **LEGEND_KWARGS)

    fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{output_prefix}.pdf", dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_calls", nargs="+", type=str, help="Input calls")
    ap.add_argument("-l", "--labels", nargs="+", type=str, help="Labels")
    ap.add_argument("-c", "--colors", nargs="+", type=str, help="Colors")
    ap.add_argument(
        "-o", "--output_prefix", default="out", type=str, help="Outfile prefix"
    )

    args = ap.parse_args()
    df_all = pl.concat(
        [
            pl.read_csv(call, separator="\t", has_header=True)
            .with_columns(lbl=pl.lit(label))
            .filter(~pl.col("name").eq("correct"))
            for call, label in zip(args.input_calls, args.labels, strict=True)
        ]
    ).with_columns(
        length=pl.col("chromEnd") - pl.col("chromStart"),
        itemRgb=pl.col("itemRgb").map_batches(rgb_to_hex, return_dtype=pl.String),
    )
    name_colors = dict(df_all.select("name", "itemRgb").iter_rows())
    df_all_grouped = df_all.group_by(["lbl", "name"]).agg(pl.col("length").sum())
    draw_bar_all(
        df_all_grouped,
        label_colors=dict(zip(args.labels, args.colors, strict=True)),
        name_colors=name_colors,
        output_prefix=args.output_prefix,
    )
    df_all_grouped.write_csv(f"{args.output_prefix}.tsv", separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
