import polars as pl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def plot(ax: Axes, df: pl.DataFrame, **kwargs) -> None:
    ax.plot(df["position"], df["first"], color="black", label="Coverage", **kwargs)
    ax.plot(df["position"], df["second"], color="red", label="Mismatch", **kwargs)
    ax.plot(
        df["position"],
        df["indel"],
        color="purple",
        label="Insertion/deletion",
        **kwargs,
    )
    for spine in ax.spines:
        if spine == "top" or spine == "right":
            ax.spines[spine].set_visible(False)


def data() -> pl.DataFrame:
    first = [9, 15, 12, 12, 11, 18, 20, 16, 11, 10]
    second = [0, 13, 1, 1, 1, 3, 5, 7, 6, 2]
    indel = [0, 2, 2, 2, 0, 1, 2, 8, 1, 2]
    return pl.DataFrame(
        {
            "position": range(len(first)),
            "first": first,
            "second": second,
            "indel": indel,
        }
    )


def data_dip() -> pl.DataFrame:
    first = [9, 15, 12, 12, 2, 18, 20, 16, 11, 10]
    second = [0, 13, 1, 1, 0, 3, 5, 7, 6, 2]
    indel = [0, 2, 2, 2, 1, 1, 2, 8, 1, 2]
    return pl.DataFrame(
        {
            "position": range(len(first)),
            "first": first,
            "second": second,
            "indel": indel,
        }
    )


def main():
    dot_kwargs = dict(marker="o", markersize=3)
    legend_kwargs = dict(bbox_to_anchor=(1, 0.5), loc="center left", frameon=False)
    df = data()
    df_filt = df.with_columns(
        pl.col("second").rolling_mean(window_size=3, center=True),
        pl.col("indel").rolling_mean(window_size=3, center=True),
    )

    # Before after filtering.
    fig, axs = plt.subplots(1, 2, layout="constrained", figsize=(10, 4.5))
    plot(axs[0], df, **dot_kwargs)
    plot(axs[1], df_filt, **dot_kwargs)

    yticks = range(0, 25, 5)
    for i, ax in enumerate(axs):
        ax: Axes
        if i == 0:
            ax.set_title("Before")
            ax.axvspan(0.5, 1.5, color="gray", alpha=0.5)
        else:
            ax.set_title("After")
        ax.set_yticks(yticks, labels=[str(tick) for tick in yticks])

    fig.supxlabel("Position")
    fig.supylabel("Sequence read depth")

    ax: Axes = axs[0]
    handles, labels = ax.get_legend_handles_labels()
    labels_handles = dict(zip(labels, handles))
    fig.legend(
        handles=labels_handles.values(), labels=labels_handles.keys(), **legend_kwargs
    )

    fig.savefig("before_after.png", bbox_inches="tight", dpi=600)

    # Dip
    fig, ax = plt.subplots(1, 1, layout="constrained", figsize=(6, 3))
    df_dip = data_dip()
    plot(ax, df_dip, **dot_kwargs)

    ax.axvspan(3.5, 4.5, color="gray", alpha=0.5)
    ax.set_title("Example: ACGACGACGT")
    ax.set_yticks(yticks, labels=[str(tick) for tick in yticks])
    ax.set_xlabel("Position")
    ax.set_ylabel("Sequence read depth")
    ax.legend(**legend_kwargs)

    fig.savefig("dip.png", bbox_inches="tight", dpi=600)

    # Peak-calling


if __name__ == "__main__":
    raise SystemExit(main())
