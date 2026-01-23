import os
import sys
import scipy.signal as sig
import polars as pl
import matplotlib.pyplot as plt

WD = os.path.dirname(__file__)


def main():
    # Homopolymer bins TSV (nt, length, count)
    bins = sys.argv[1]
    df = pl.read_csv(bins, separator="\t", new_columns=["nt", "len", "count"])
    df = df.group_by(["len"]).agg(pl.col("count").sum()).sort("len")

    fig, ax = plt.subplots()

    peaks = sig.find_peaks(df["count"], prominence=5, height=1000)

    ax.bar(df["len"], df["count"])
    for peak, ht in zip(df["len"][peaks[0]], peaks[1]["peak_heights"]):
        ax.axvline(peak, linestyle="dotted", color="black")
        ax.annotate(f"{peak}", xy=(peak, ht))

    fig.savefig(os.path.join(WD, "out.png"))


if __name__ == "__main__":
    raise SystemExit(main())
