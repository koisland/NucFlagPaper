import argparse
from matplotlib.colors import rgb2hex
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes
from matplotlib.text import Text


# https://stackoverflow.com/a/67594395
def update_bars(ax: Axes, round_to: int):
    for c in ax.containers:
        labels = [f"{round(v.get_height() / 1_000_000, round_to)}" for v in c]
        ax.bar_label(c, labels=labels, label_type="edge")
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
    colors: dict[str, str],
    name_colors: dict[str, str],
    col_order: list[str],
    outfile: str,
):
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(col_order),
        figsize=(16, 4),
        layout="constrained",
        sharex=True,
        sharey=True,
    )
    max_length = df_all["length"].max()
    yticks = list(range(0, int(max_length), 1_000_000))
    yticklabels = [f"{tick / 1_000_000}" for tick in yticks]
    for i, col in enumerate(col_order):
        ax: Axes = axes[i]
        df_col = df_all.filter(pl.col("lbl") == col)

        sns.barplot(
            data=df_col,
            x="name",
            y="length",
            hue="lbl",
            legend=None,
            palette=colors,
            ax=ax,
        )

        update_bars(ax, round_to=2)
        for lbl in ax.get_xticklabels():
            lbl: Text
            # Color and stroke around text
            lbl.set_color(name_colors[lbl.get_text()])
            lbl.set_path_effects([pe.Stroke(linewidth=0.2, foreground="black")])

        ax.set_yticks(yticks, yticklabels)
        ax.set_title(col)
        ax.tick_params(axis="x", rotation=45)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.supylabel("Length (Mbp)")
    fig.savefig(outfile, dpi=600, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_calls", nargs="+", type=str, help="Input calls")
    ap.add_argument("-l", "--labels", nargs="+", type=str, help="Labels")
    ap.add_argument("-c", "--colors", nargs="+", type=str, help="Colors")
    ap.add_argument("-o", "--outfile", default="out.png", type=str, help="Outfile")

    args = ap.parse_args()

    colors = dict(zip(args.labels, args.colors, strict=True))

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
    draw_bar_all(
        df_all.group_by(["lbl", "name"]).agg(pl.col("length").sum()),
        colors=colors,
        col_order=args.labels,
        name_colors=name_colors,
        outfile=args.outfile,
    )


if __name__ == "__main__":
    raise SystemExit(main())
