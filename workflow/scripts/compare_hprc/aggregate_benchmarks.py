import os
import glob
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes


def read_files(fglob: str) -> pl.DataFrame:
    dfs: list[pl.DataFrame] = []
    for file in glob.glob(fglob):
        _, ext = os.path.splitext(file)

        fname = os.path.basename(file).replace(ext, "")
        _, sm, release = fname.rsplit("_", 2)
        release = release.replace("R", "Release ")
        dfs.append(
            pl.read_csv(
                file,
                separator="\t",
                has_header=True,
            ).with_columns(
                sm=pl.lit(sm),
                release=pl.lit(release),
            )
        )
    return pl.concat(dfs).with_columns(
        **{
            "Wall Time (hours)": pl.col("s") / 3600,
            "CPU Time (hours)": pl.col("cpu_time") / 3600,
            "Max Resident Set Size (GB)": pl.col("max_rss") / 1024,
        }
    )


def plot_agg(df: pl.DataFrame, output_prefix: str, metrics: list[str]):
    for metric in metrics:
        metric_fs = metric.replace(" ", "_").replace("(", "_").replace(")", "")
        fig, ax = plt.subplots()
        ax: Axes
        sns.violinplot(
            data=df,
            x="release",
            y=metric,
            hue="release",
            hue_order=["Release 1", "Release 2"],
            palette={"Release 1": "red", "Release 2": "blue"},
            inner="quart",
            ax=ax,
            alpha=0.6,
        )
        sns.stripplot(
            data=df,
            x="release",
            y=metric,
            hue="release",
            edgecolor="black",
            linewidth=0.5,
            hue_order=["Release 1", "Release 2"],
            palette={"Release 1": "red", "Release 2": "blue"},
            ax=ax,
        )
        ax.set_xlabel(None)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        for line in ax.lines:
            data = line.get_data()
            ax.text(
                data[0][data[0].nonzero()][0],
                data[1][0],
                f"{data[1][0]:.1f}",
                size=10,
                zorder=2,
                path_effects=[pe.withStroke(linewidth=2, foreground="w")],
            )

        prefix = f"{output_prefix}_{metric_fs}"
        fig.savefig(f"{prefix}.png", bbox_inches="tight", dpi=300)
        fig.savefig(f"{prefix}.pdf", bbox_inches="tight", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--bmk_dir", default="benchmarks/hprc")
    ap.add_argument(
        "-o", "--output_prefix", default="results/hprc/benchmarks/hprc_runtime_metrics"
    )
    args = ap.parse_args()
    dname = os.path.dirname(args.output_prefix)
    os.makedirs(dname, exist_ok=True)

    df = read_files(os.path.join(args.bmk_dir, "run_nucflag*.tsv"))
    plot_agg(
        df,
        args.output_prefix,
        ["Wall Time (hours)", "CPU Time (hours)", "Max Resident Set Size (GB)"],
    )
    df.write_csv(f"{args.output_prefix}.tsv", separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
