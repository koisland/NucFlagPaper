import argparse
import polars as pl
import seaborn as sns

from matplotlib.axes import Axes


# https://stackoverflow.com/a/67594395
def update_bars(ax: Axes, round_to: int):
    for c in ax.containers:
        labels = [f"{round(v.get_height() / 1_000_000, round_to)}" for v in c]
        ax.bar_label(c, labels=labels, label_type="edge")

    ylim = ax.get_ylim()
    yticks = range(int(ylim[0]), int(ylim[1]), 10_000_000)
    ax.set_yticks(yticks, [f"{round(tick / 1_000_000, round_to)}" for tick in yticks])
    ax.set_xlabel(None)


def draw_bar_all(df_all: pl.DataFrame, colors: dict[str, str], outfile: str):
    g = sns.catplot(
        data=df_all,
        x="name",
        y="length",
        col="lbl",
        hue="lbl",
        sharex=False,
        kind="bar",
        legend=None,
        aspect=1.5,
        palette=colors,
    )
    for ax in g.axes.ravel():
        update_bars(ax, round_to=2)

    g.set_titles(row_template="{row_name}", col_template="{col_name}")

    g.tick_params(axis="x", rotation=45)
    g.set_ylabels("Length (Mbp)")
    g.savefig(outfile, dpi=600, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_calls", nargs="+", type=str, help="Input calls")
    ap.add_argument("-l", "--labels", nargs="+", type=str, help="Labels")
    ap.add_argument("-c", "--colors", nargs="+", type=str, help="Colors")
    ap.add_argument("-o", "--outfile", default="out.png", type=str, help="Outfile")

    args = ap.parse_args()

    colors = dict(zip(args.labels, args.colors, strict=True))

    df_all = (
        pl.concat(
            [
                pl.read_csv(call, separator="\t", has_header=True).with_columns(
                    lbl=pl.lit(label)
                )
                for call, label in zip(args.calls, args.labels, strict=True)
            ]
        )
        .with_columns(length=pl.col("chromEnd") - pl.col("chromStart"))
        .group_by(["lbl", "name"])
        .agg(pl.col("length").sum())
    )

    draw_bar_all(df_all, colors=colors, outfile=args.outfile)


if __name__ == "__main__":
    raise SystemExit(main())
