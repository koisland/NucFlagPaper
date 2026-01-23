import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    title=None,
    fancybox=False,
    frameon=False,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_qvs", nargs="+", required=True)
    ap.add_argument("-r", "--reported_qvs", required=True)
    ap.add_argument("-l", "--labels", nargs="+", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    dfs_qv: list[pl.DataFrame] = []
    for bed, label in zip(args.input_qvs, args.labels, strict=True):
        # Assumes chromosomes in order
        df = (
            pl.read_csv(bed, separator="\t")
            .with_columns(label=pl.lit(label), chrom=(pl.col("#chrom").rle_id() + 1))
            .filter(pl.col("chrom") < 6)
            .select("chrom", "label", "QV")
        )
        dfs_qv.append(df)

    df_qv = (
        pl.concat(dfs_qv)
        .join(
            pl.read_csv(args.reported_qvs, has_header=True, separator="\t"),
            how="left",
            on=["chrom", "label"],
        )
        .rename({"QV": "QV_NucFlag", "QV_right": "QV_Merqury"})
        .unpivot(
            on=["QV_NucFlag", "QV_Merqury"],
            index=["chrom", "label"],
            variable_name="Method",
            value_name="QV",
        )
    )
    fig, axes = plt.subplots(
        nrows=1,
        ncols=df_qv["chrom"].n_unique(),
        figsize=(16, 4),
        layout="constrained",
        sharex=True,
        sharey=True,
    )
    for i, chrom in enumerate(sorted(df_qv["chrom"].unique())):
        ax: Axes = axes[i]
        df_chrom = df_qv.filter(pl.col("chrom") == chrom)

        sns.barplot(
            data=df_chrom,
            x="Method",
            y="QV",
            hue="label",
            order=["QV_Merqury", "QV_NucFlag"],
            palette="colorblind",
            ax=ax,
            legend="full",
        )
        for c in ax.containers:
            labels = [f"{round(v.get_height(), 1)}" for v in c]
            ax.bar_label(c, labels=labels, label_type="edge")

        ax.set_xlabel(None)
        ax.set_title(f"Chromosome {chrom}")

        if i != 0:
            ax.get_legend().remove()

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    ax_legend: Axes = axes[0]
    handles, labels = ax_legend.get_legend_handles_labels()
    fig.legend(
        labels=labels,
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncols=len(labels),
        **LEGEND_KWARGS,
    )
    ax_legend.get_legend().remove()
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
