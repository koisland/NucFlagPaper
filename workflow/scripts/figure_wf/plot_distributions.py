import os
from matplotlib.patches import Patch
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def main():
    typs = ["cov", "indel", "mapq", "mismatch", "softclip"]
    fig, axes = plt.subplots(
        nrows=1, ncols=len(typs), layout="constrained", figsize=(12, 3)
    )

    for i, typ in enumerate(typs):
        ax: Axes = axes[i]

        df = pl.read_csv(
            os.path.join(
                "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/HG002_ont_chr10_5Mbp",
                f"chr10_PATERNAL_{typ}.bw.bg",
            ),
            has_header=False,
            separator="\t",
            new_columns=["chrom", "st", "end", "value"],
        )
        if not (typ == "cov" or typ == "mapq"):
            df = df.filter(pl.col("value") != 0.0)

        min_value = df["value"].min()
        max_value = df["value"].max()
        assert isinstance(min_value, float)
        assert isinstance(max_value, float)

        median = df["value"].median()
        adj_stdev = (df["value"] - median).abs().median()
        assert isinstance(adj_stdev, float)
        assert isinstance(median, float)
        if adj_stdev == 0.0:
            adj_stdev = (df["value"] - df["value"].mean()).abs().mean()
            assert isinstance(adj_stdev, float)
        else:
            adj_stdev = adj_stdev * 1.486

        print(typ, median, adj_stdev)
        ax.hist(df["value"], bins=100)
        ax.set_xlabel(typ)
        ax.set_ylabel("Count")
        ax.axvline(x=median, color="black", linestyle="dotted")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        for i in range(1, 4):
            x_iminus_zscore = median - (adj_stdev * i)
            x_iplus_zscore = median + (adj_stdev * i)
            if not x_iminus_zscore < min_value:
                ax.axvline(x=x_iminus_zscore, color="red", linestyle="dotted")
            if not x_iplus_zscore > max_value:
                ax.axvline(x=x_iplus_zscore, color="red", linestyle="dotted")

        if typ == "cov" or typ == "mapq":
            continue

        ax.set_yscale("log")

    fig.legend(
        handles=[
            Patch(label="Median", color="black"),
            Patch(label="Adjusted\nz-score", color="red"),
        ],
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
    )
    fig.savefig("dist.png", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
