import sys
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

RELEASES = (1, 2)
RELEASES_COLORS = ("red", "blue")
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)


def main():
    r1 = sys.argv[1]
    r2 = sys.argv[2]
    plot = sys.argv[3]

    r1_fmt_repl_map = {"X": "", "-": "0"}
    df_r1 = (
        pl.read_csv(r1, separator="\t", has_header=True)
        .select("Sample", "ONT >100kb")
        .with_columns(
            pl.col("ONT >100kb").str.replace_many(r1_fmt_repl_map).cast(pl.UInt64),
            release=pl.lit(1),
        )
    )
    samples_r1 = df_r1["Sample"]
    df_r2 = (
        pl.read_csv(r2, separator=",", has_header=True)
        .filter(pl.col("sample_ID").is_in(samples_r1))
        .group_by(["sample_ID"])
        .agg(
            **{
                "ONT >100kb": pl.col("100kb+").sum().cast(pl.UInt64),
            }
        )
        .rename({"sample_ID": "Sample"})
        .with_columns(release=pl.lit(2))
    )
    df_both = pl.concat([df_r1, df_r2])
    df_wide_both = df_both.pivot(
        on="release", index="Sample", values="ONT >100kb"
    ).with_columns(ratio=pl.col("2") / pl.col("1"))
    print(f"Median increase: {df_wide_both['ratio'].median()}", file=sys.stderr)
    fig, ax = plt.subplots(figsize=(16, 8), layout="constrained")
    sns.barplot(
        df_both,
        x="Sample",
        y="ONT >100kb",
        hue="release",
        palette=dict(zip(RELEASES, RELEASES_COLORS)),
        ax=ax,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_xlabel(None)
    ax.set_ylabel("Coverage of ONT reads >100kb (X)")
    sns.move_legend(ax, loc="upper right", title="HPRC release", **LEGEND_KWARGS)
    ax.tick_params(axis="x", labelrotation=45.0)
    fig.savefig(plot, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
