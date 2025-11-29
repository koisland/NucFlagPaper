import os
import re
import argparse
import polars as pl
import seaborn as sns

RGX_ALIGNER_FMT = re.compile(r"nucflag_HG002_(hifi|ont_r10)_(.*?)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infiles",
        nargs="+",
        metavar="{method}_{sample}_{dtype}",
        help="Input summary files. Expects header.",
    )
    ap.add_argument("-o", "--output_dir", default=".", help="Output dir.")

    args = ap.parse_args()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    dfs = []
    for file in args.infiles:
        dtype, aligner = RGX_ALIGNER_FMT.search(
            os.path.splitext(os.path.basename(file))[0]
        ).groups()
        dfs.append(
            pl.read_csv(
                file,
                separator="\t",
                schema_overrides={"downsample": pl.Float32},
                null_values=["None"],
            ).with_columns(method=pl.lit(aligner))
        )

    df = (
        pl.concat(dfs)
        .sort("method", "dtype", "mtype")
        .with_columns(
            downsample=pl.when(pl.col("downsample").is_null())
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
            palette="Paired",
            height=height,
            aspect=aspect_ratio,
        )
        for ax in g.axes.flat:
            ax.set_ylim(0, 1.0)

        g.set_axis_labels(label, "Precision (%)")
        g.savefig(f"{output_dir}/precision_{var}.png")

        g = sns.catplot(
            data=df,
            kind="bar",
            x=var,
            y="recall",
            row="dtype",
            col="mtype",
            hue="method",
            errorbar="sd",
            palette="Paired",
            height=height,
            aspect=aspect_ratio,
        )
        for ax in g.axes.flat:
            ax.set_ylim(0, 1.0)

        g.set_axis_labels(label, "Recall (%)")
        g.savefig(f"{output_dir}/recall_{var}.png")


if __name__ == "__main__":
    raise SystemExit(main())
