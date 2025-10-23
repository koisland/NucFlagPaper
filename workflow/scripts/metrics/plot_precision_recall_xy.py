import os
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.axes import Axes

DF_READ_DEPTH = pl.DataFrame(
    {
        "sample": ["HG002", "HG002", "A_thaliana"],
        "coverage": [49.9671, 21.1335, 148.523],
        "dtype": ["hifi", "ont_r10", "hifi"],
    }
)
TOOL_VERSIONS = {
    "nucflag": "NucFlag v1.0.0",
    "inspector": "Inspector v1.3",
    "flagger": "HMM-Flagger v1.1.0",
}
TOOL_COLORS = {
    "NucFlag v1.0.0": "purple",
    "Inspector v1.3": "teal",
    "HMM-Flagger v1.1.0": "magenta",
}
MISASSEMBLY_NAMES = {
    "misjoin": "Misjoin",
    "false_duplication": "False Duplication",
    "inversion": "Inversion",
}
DTYPE_NAMES = {"hifi": "PacBio HiFi", "ont_r10": "ONT (R10)"}


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

    dfs: list[pl.DataFrame] = []
    for file in args.infiles:
        method, sm, dtype = os.path.splitext(os.path.basename(file))[0].split("_", 2)
        dfs.append(
            pl.read_csv(
                file,
                separator="\t",
                schema_overrides={"downsample": pl.Float32},
                null_values=["None"],
            )
            .with_columns(Method=pl.lit(TOOL_VERSIONS[method]))
            .rename({"mtype": "Misassembly"})
        )

    df = (
        pl.concat(dfs)
        .sort("Method", "dtype", "Misassembly")
        .with_columns(
            perc_coverage=pl.when(pl.col("downsample").is_null())
            .then(pl.lit("1.0"))
            .otherwise(pl.col("downsample"))
            .cast(pl.Float32)
        )
        .join(DF_READ_DEPTH, on=["sample", "dtype"])
        .with_columns(coverage=(pl.col("perc_coverage") * pl.col("coverage")).round(1))
        .with_columns(pl.col("recall") * 100, pl.col("precision") * 100)
    )
    remove_spines = ["top", "right"]
    dtype_order = ["hifi", "ont_r10"]
    misassembly_order = ["misjoin", "inversion", "false_duplication"]
    for var, xlabel in (("len", "Length (bp)"), ("coverage", "Read depth (X)")):
        dfs_dtype = df.sort(var).partition_by(["dtype", "Misassembly"], as_dict=True)

        for metric in ("precision", "recall"):
            fig, axes = plt.subplots(
                nrows=2,
                ncols=3,
                sharex=var == "len",
                layout="constrained",
                figsize=(16, 8),
                height_ratios=[0.5, 0.5],
            )
            labels_handles = {}
            facet_order = [
                (row, col, dtype, misassembly)
                for row, dtype in enumerate(dtype_order)
                for col, misassembly in enumerate(misassembly_order)
            ]

            fig.supxlabel(xlabel)

            for row, col, dtype, misassembly in facet_order:
                df_dtype = dfs_dtype.get((dtype, misassembly))
                if not isinstance(df_dtype, pl.DataFrame):
                    continue

                ax: Axes = axes[row, col]

                ax.set_ylim(0.0, 100.0)
                ax.set_xlim(df_dtype[var].min(), df_dtype[var].max())

                if var == "len":
                    sns.lineplot(
                        data=df_dtype,
                        x=var,
                        y=metric,
                        hue="Method",
                        palette=TOOL_COLORS,
                        style="Misassembly",
                        markers=True,
                        dashes=False,
                        ax=ax,
                    )
                    ax.set_xscale("log")
                    values = df_dtype[var].unique().sort().to_list()
                    ax.set_xticks(values, [str(val) for val in values])
                else:
                    sns.barplot(
                        data=df_dtype,
                        x=var,
                        y=metric,
                        hue="Method",
                        palette=TOOL_COLORS,
                        ax=ax,
                        # Reverse order.
                        hue_order=TOOL_COLORS.keys(),
                        order=df_dtype[var].unique().sort(descending=True),
                    )

                handles, labels = ax.get_legend_handles_labels()
                for handle, label in zip(handles, labels):
                    if label in TOOL_VERSIONS.values():
                        labels_handles[label] = handle

                ax.legend().set_visible(False)

                for spine in remove_spines:
                    ax.spines[spine].set_visible(False)

                ax.set_xlabel(None)
                ax.set_ylabel(None)

                if col == 0:
                    ax.set_ylabel(
                        f"{metric.capitalize()} (%)\n({DTYPE_NAMES[dtype]})",
                        rotation=0,
                        ha="right",
                        ma="center",
                        fontsize="large",
                    )
                if row == 0:
                    ax.set_title(MISASSEMBLY_NAMES[misassembly])

            fig.legend(
                labels_handles.values(),
                labels_handles.keys(),
                loc="outside center right",
            )
            fig.savefig(f"{output_prefix}{metric}_{var}.png")
            fig.clear()


if __name__ == "__main__":
    raise SystemExit(main())
