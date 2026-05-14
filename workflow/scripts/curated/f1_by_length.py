import sys
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

plt.rcParams["font.family"] = "Arial"


def main():
    infile = sys.argv[1]
    outprefix = sys.argv[2]

    df = (
        pl.read_csv(infile, separator="\t", has_header=True)
        .with_columns(length=(pl.col("end") - 5) - (pl.col("st") + 5))
        # .with_columns(
        #     group=pl.when(pl.col("length").is_between(1, 5, closed="left"))
        #     .then(pl.lit("1-5"))
        #     .when(pl.col("length").is_between(5, 10, closed="left"))
        #     .then(pl.lit("5-10"))
        #     .when(pl.col("length").is_between(10, 50, closed="left"))
        #     .then(pl.lit("10-50"))
        #     .when(pl.col("length").is_between(50, 100, closed="left"))
        #     .then(pl.lit("50-100"))
        #     .when(pl.col("length").is_between(100, 500, closed="left"))
        #     .then(pl.lit("100-500"))
        #     .when(pl.col("length").is_between(500, 1000, closed="left"))
        #     .then(pl.lit("500-1000"))
        #     .when(pl.col("length").is_between(1000, 5000, closed="left"))
        #     .then(pl.lit("1000-5000"))
        #     .when(pl.col("length").is_between(5000, 10000, closed="left"))
        #     .then(pl.lit("5000-10000"))
        #     .otherwise(pl.lit(">10000"))
        # )
        .group_by(["length"])
        .agg(pl.col("type").value_counts())
        .explode("type")
        .unnest("type")
        .pivot(index="length", values="count", on="type")
        .fill_null(pl.lit(0))
        .with_columns(
            recall=pl.col("true_positive")
            / (pl.col("true_positive") + pl.col("false_negative")),
            precision=pl.col("true_positive")
            / (pl.col("true_positive") + pl.col("false_positive")),
        )
        .sort(by="length")
        .with_columns(
            f1=2
            * (
                (pl.col("recall") * pl.col("precision"))
                / (pl.col("recall") + pl.col("precision"))
            )
        )
        .fill_nan(0)
    )
    fig, ax = plt.subplots()
    ax: Axes
    # xticks, xlabels = zip(*df.select("group_st", "group").iter_rows())
    sns.pointplot(data=df, x="length", y="f1", ax=ax, label="F1")
    # sns.pointplot(
    #     data=df,
    #     x="length",
    #     y="recall",
    #     ax=ax,
    #     color="orange",
    #     label="recall"
    # )
    # sns.pointplot(
    #     data=df,
    #     x="length",
    #     y="precision",
    #     ax=ax,
    #     color="green",
    #     label="precision"
    # )
    ax.set_ylim(0, 1)
    ax.legend()
    # ax.set_xticks(xticks, xlabels)
    fig.savefig(f"{outprefix}.png")


if __name__ == "__main__":
    raise SystemExit(main())
