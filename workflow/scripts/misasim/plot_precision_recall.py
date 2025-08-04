import os
import argparse
import polars as pl
import seaborn as sns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infiles",
        nargs="+",
        metavar="{method}_{sample}_{dtype}",
        help="Input summary files. Expects header.",
    )
    ap.add_argument("-o", "--output_prefix", default="./", help="Output prefix.")

    args = ap.parse_args()
    output_prefix = args.output_prefix
    output_dir = os.path.dirname(output_prefix)
    os.makedirs(output_dir, exist_ok=True)

    dfs = []
    for file in args.infiles:
        method, sm, dtype = os.path.splitext(os.path.basename(file))[0].split("_", 2)
        dfs.append(
            pl.read_csv(file, separator="\t").with_columns(method=pl.lit(method))
        )

    df = (
        pl.concat(dfs)
        .sort("method", "dtype", "mtype")
        .with_columns(
            downsample=pl.when(pl.col("downsample") == "None")
            .then(pl.lit("0.0"))
            .otherwise(pl.col("downsample"))
            .cast(pl.Float32)
        )
        .with_columns(pl.col("downsample") * 100)
    )

    height = 3
    aspect_ratio = 1.5
    for var, label in (("len", "Length (bp)"), ("downsample", "Downsample ratio (%)")):
        df = df.sort(var)

        g = sns.catplot(
            data=df,
            kind="bar",
            x=var,
            y="precision",
            row="dtype",
            col="mtype",
            hue="method",
            errorbar="sd",
            palette="pastel",
            height=height,
            aspect=aspect_ratio,
        )
        for ax in g.axes.flat:
            ax.set_ylim(0, 1.0)

        g.set_axis_labels(label, "Precision (%)")
        g.savefig(f"{output_prefix}precision_{var}.png")

        g = sns.catplot(
            data=df,
            kind="bar",
            x=var,
            y="recall",
            row="dtype",
            col="mtype",
            hue="method",
            errorbar="sd",
            palette="pastel",
            height=height,
            aspect=aspect_ratio,
        )
        for ax in g.axes.flat:
            ax.set_ylim(0, 1.0)

        g.set_axis_labels(label, "Recall (%)")
        g.savefig(f"{output_prefix}recall_{var}.png")


if __name__ == "__main__":
    raise SystemExit(main())
