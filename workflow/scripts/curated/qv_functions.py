import sys
import matplotlib.pyplot as plt

from math import log10
from matplotlib.axes import Axes


# Merqury's QV metric.
# https://github.com/marbl/merqury/blob/1ad7c328a0098ba7a5edb85a5f6a18fe62869e46/eval/qv.sh#L82
# $ echo "1 1670" | awk -v k=31 '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'
def qv_merqury(
    kmer_asm: int, kmer_total: int, kmer_size: int = 31, ndigits: int = 3
) -> float:
    return round(
        -10 * log10(1 - (1 - kmer_asm / kmer_total) ** (1 / kmer_size)) / log10(10),
        ndigits,
    )


# Based on Inspector's QV metric.
# https://github.com/ChongLab/Inspector/blob/0e08f882181cc0e0e0fa749cd87fb74a278ea0f0/inspector.py#L184
# $ python -c "import math,sys; e,t = sys.argv[1:3]; print(-10 * math.log10(int(e) / int(t)))" 1 1670
def qv_nucflag(bp_err: int, bp_total: int, ndigits: int = 3) -> float:
    return round(-10 * log10(bp_err / bp_total), ndigits)


def main():
    output = sys.argv[1]

    """
    +----------------------------+
    |            |merqury|nucflag|
    +------------+-------+-------+
    |bp_err      |       |       |
    +----------------------------+
    """
    num_increment = 8
    kmer_size = 31

    sizes = [10**i for i in reversed(range(num_increment + 1))]
    max_size = sizes[0]
    total_kmers = max_size - kmer_size + 1
    kmers_sizes = [min(size, total_kmers) for size in sizes]

    all_qv_nucflag_bp_err = [(bp, abs(qv_nucflag(bp, max_size))) for bp in sizes]
    all_qv_merqury_bp_err = [
        (asm_kmers, abs(qv_merqury(asm_kmers, total_kmers, kmer_size)))
        for asm_kmers in kmers_sizes
    ]
    fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, layout="constrained")
    fig: plt.Figure
    ax: Axes
    axes_info = [
        (
            axes[0],
            all_qv_merqury_bp_err,
            "Merqury",
            "Number of assembly only kmers (k=31)",
            "blue",
        ),
        (axes[1], all_qv_nucflag_bp_err, "NucFlag", "Length of calls (bp)", "red"),
    ]
    for ax, data, title, xlabel, color in axes_info:
        ax.set_xscale("log")
        ax.set_title(title)
        ax.set_xlabel(xlabel)

        for val, qv in data:
            ax.scatter(val, qv, color=color)
            ax.annotate(str(qv), (val, qv), xytext=(0, 10), textcoords="offset points")

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.supylabel("QV")
    fig.savefig(output)


if __name__ == "__main__":
    raise SystemExit(main())
