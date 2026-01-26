import os
import glob
import argparse
import matplotlib.colors
import polars as pl
import seaborn as sns

RGX_DOWNSAMPLE = r"(?<downsample>0\.33|0\.5)"
RGX_STEP_SM_DTYPE = r"(?<step>^.*)_(?<sample>.*?)_(?<dtype>hifi|ont_r10|ont_r9)"
RGX_MTYPE_NUM_LEN = (
    r"(?<mtype>false_duplication|misjoin|inversion)-(?<num>\d+)-(?<len>\d+)"
)
TOOL_NAMES = [
    "NucFlag v1.0\n(CRAM)",
    "Inspector v1.3",
    "HMM-Flagger v1.1.0",
]
TOOL_COLORS = ["purple", "teal", "magenta"]


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
            pl.col("max_rss").max() / 1024,
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
    height = 3
    aspect_ratio = 1.5
    colors = dict(
        zip(
            TOOL_NAMES,
            [matplotlib.colors.to_rgba(color, alpha=0.5) for color in TOOL_COLORS],
        )
    )
    g = sns.catplot(
        data=df,
        x="Method",
        y=y,
        row="Data",
        col="Misassembly",
        order=order,
        hue_order=order,
        palette=colors,
        kind="box",
        hue="Method",
        height=height,
        aspect=aspect_ratio,
    )
    # https://stackoverflow.com/a/69398767
    g.map(
        sns.stripplot,
        "Method",
        y,
        "Method",
        order=order,
        hue_order=order,
        palette=colors,
        size=3,
        linewidth=1,
    )

    g.set(xlabel=None)
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nucflag_bmks_dir", default="benchmarks/misasim/nucflag")
    ap.add_argument(
        "-i", "--inspector_bmks_dir", default="benchmarks/misasim/inspector"
    )
    ap.add_argument("-f", "--flagger_bmks_dir", default="benchmarks/misasim/flagger")
    ap.add_argument("-o", "--output_dir", default="results/misasim/benchmark_plots")
    args = ap.parse_args()

    df = (
        pl.concat(
            [
                read_files(os.path.join(args.flagger_bmks_dir, "*.txt")).with_columns(
                    method=pl.lit("HMM-Flagger v1.1.0"),
                    step=pl.lit("run_flagger"),
                ),
                read_files(os.path.join(args.inspector_bmks_dir, "*.txt")).with_columns(
                    method=pl.lit("Inspector v1.3")
                ),
                read_files(
                    os.path.join(args.nucflag_bmks_dir, "run_nucflag*.tsv")
                ).with_columns(method=pl.lit("NucFlag v1.0\n(CRAM)")),
            ]
        )
        .sort("method", "dtype", "mtype")
        .rename(
            {
                "dtype": "Data",
                "mtype": "Misassembly",
                "method": "Method",
                "max_rss": "Maximum RSS (GB)",
                "hours": "Wall Time (hours)",
            }
        )
    )

    os.makedirs(args.output_dir, exist_ok=True)

    df_summary = df.group_by(["Method", "Data"]).agg(
        pl.col("Wall Time (hours)").median(), pl.col("Maximum RSS (GB)").median()
    )
    df_summary.write_csv(os.path.join(args.output_dir, "summary.csv"))

    g_rss = plot_agg(df.to_pandas(), "Maximum RSS (GB)", TOOL_NAMES)
    g_rss.savefig(os.path.join(args.output_dir, "rss_gb.png"), dpi=600)
    g_time = plot_agg(df.to_pandas(), "Wall Time (hours)", TOOL_NAMES)
    g_time.savefig(os.path.join(args.output_dir, "hours.png"), dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
