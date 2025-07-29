import os
import glob
import argparse
import polars as pl
import seaborn as sns

RGX_DOWNSAMPLE = r"(?<downsample>0\.33|0\.5)"
RGX_STEP_SM_DTYPE = r"(?<step>^.*)_(?<sample>.*?)_(?<dtype>hifi|ont_r10|ont_r9)"
RGX_MTYPE_NUM_LEN = (
    r"(?<mtype>false_duplication|misjoin|inversion)-(?<num>\d+)-(?<len>\d+)"
)


def read_files(fglob: str) -> pl.DataFrame:
    dfs: list[pl.DataFrame] = []
    for file in glob.glob(fglob):
        _, ext = os.path.splitext(file)

        fname = os.path.basename(file).replace(ext, "")
        dfs.append(
            pl.read_csv(
                file,
                separator="\t",
                has_header=True,
            ).with_columns(fname=pl.lit(fname))
        )
    df = (
        pl.concat(dfs)
        .with_columns(
            downsample=pl.col("fname").str.extract_groups(RGX_DOWNSAMPLE),
            step_sm_dtype=pl.col("fname").str.extract_groups(RGX_STEP_SM_DTYPE),
            mtype_num_len=pl.col("fname").str.extract_groups(RGX_MTYPE_NUM_LEN),
        )
        .unnest("downsample", "step_sm_dtype", "mtype_num_len")
        .group_by(["sample", "dtype", "downsample", "mtype", "num", "len"])
        .agg(
            pl.col("step").first(),
            pl.col("max_rss").max(),
            pl.col("max_vms").max(),
            pl.col("max_uss").max(),
            pl.col("max_pss").max(),
            pl.col("io_in").mean(),
            pl.col("io_out").mean(),
            pl.col("mean_load").mean(),
            pl.col("cpu_time").sum(),
            minutes=pl.col("s").sum() / 60,
            hours=pl.col("s").sum() / 3600,
        )
    )
    return df


def plot_agg(df: pl.DataFrame, y: str, order: list[str]) -> sns.FacetGrid:
    g = sns.catplot(
        data=df,
        x="method",
        y=y,
        row="dtype",
        col="mtype",
        order=order,
        hue_order=order,
        kind="box",
        hue="method",
    )
    # https://stackoverflow.com/a/69398767
    g.map(
        sns.swarmplot,
        "method",
        y,
        "method",
        order=order,
        hue_order=order,
        palette=sns.color_palette()[0:len(order)],
        alpha=0.6,
        size=3,
        linewidth=1
    )
    return g

def main():
    ap = argparse.ArgumentParser()
    # /project/logsdon_shared/projects/Keith/NucFlagPaper/benchmarks/misasim/nucflag
    ap.add_argument("-n", "--nucflag_bmks_dir", default="benchmarks/misasim/nucflag")
    # /project/logsdon_shared/projects/Keith/NucFlagPaper/benchmarks/misasim/inspector
    ap.add_argument("-i", "--inspector_bmks_dir", default="benchmarks/misasim/inspector")
    # /project/logsdon_shared/projects/Keith/NucFlagPaper/benchmarks/misasim/flagger
    ap.add_argument("-f", "--flagger_bmks_dir", default="benchmarks/misasim/flagger")
    ap.add_argument("-o", "--output_dir", default="results/benchmarks")
    args = ap.parse_args()

    df = pl.concat(
        [
            read_files(os.path.join(args.flagger_bmks_dir, "*.txt")).with_columns(
                method=pl.lit("flagger"),
                step=pl.lit("run_flagger"),
            ),
            read_files(os.path.join(args.inspector_bmks_dir, "*.txt")).with_columns(
                method=pl.lit("inspector")
            ),
            read_files(
                os.path.join(args.nucflag_bmks_dir, "run_nucflag*.tsv")
            ).with_columns(method=pl.lit("nucflag")),
        ]
    ).to_pandas()

    order = ["flagger", "inspector", "nucflag"]
    g_rss = plot_agg(df, "max_rss", order)
    g_rss.savefig(os.path.join(args.output_dir, "rss.png"))
    g_time = plot_agg(df, "minutes", order)
    g_time.savefig(os.path.join(args.output_dir, "minutes.png"))


if __name__ == "__main__":
    raise SystemExit(main())
