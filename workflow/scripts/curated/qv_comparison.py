import argparse
import scipy.stats

from matplotlib.axes import Axes
from matplotlib.lines import Line2D

import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Source - https://stackoverflow.com/a
# Posted by ImportanceOfBeingErnest, modified by community. See post 'Timeline' for change history
# Retrieved 2025-12-24, License - CC BY-SA 4.0
class SeabornFig2Grid:
    def __init__(
        self, seaborngrid, fig: plt.Figure, subplot_spec: gridspec.SubplotSpec
    ):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or isinstance(
            self.sg, sns.axisgrid.PairGrid
        ):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """Move PairGrid or Facetgrid"""
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n, m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i, j], self.subgrid[i, j])

    def _movejointgrid(self):
        """Move Jointgrid"""
        h = self.sg.ax_joint.get_position().height
        h2 = self.sg.ax_marg_x.get_position().height
        r = int(np.round(h / h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(
            r + 1, r + 1, subplot_spec=self.subplot
        )

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax: Axes, gs: gridspec.SubplotSpec):
        # https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure = self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-l",
        "--labels",
        nargs="+",
        required=True,
        help="Labels of equal number of x and y.",
        type=str,
    )
    ap.add_argument(
        "-c",
        "--colors",
        nargs="+",
        required=True,
        help="Colors of equal number of x and y.",
        type=str,
    )
    ap.add_argument(
        "-x",
        "--qv_x",
        nargs="+",
        required=True,
        help="QV of x.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-y",
        "--qv_y",
        nargs="+",
        required=True,
        help="QV of y.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument("--xlabel", default="QV (Merqury)", help="x-label", type=str)
    ap.add_argument("--ylabel", default="QV (NucFlag)", help="y-label", type=str)
    ap.add_argument(
        "--omit-acrocentrics",
        dest="omit_acrocentrics",
        action="store_true",
        help="Omit acrocentric chromosomes.",
    )
    ap.add_argument("-o", "--output", required=True, help="Output plot.", type=str)

    # x - nucflag, y - merqury
    args = ap.parse_args()

    def annotate(ax, data, **kws):
        r, p = scipy.stats.pearsonr(data["qv_x"], data["qv_y"])
        ax.text(0.05, 0.8, "r={:.2f}, p={:.2g}".format(r, p), transform=ax.transAxes)

    figs = []
    for label, color, qv_x, qv_y in zip(
        args.labels, args.colors, args.qv_x, args.qv_y, strict=True
    ):
        # Expects chrom, qv
        df_x = pl.read_csv(
            qv_x,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_x"],
            comment_prefix="#",
            has_header=False,
        )
        df_y = pl.read_csv(
            qv_y,
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "qv_y"],
            comment_prefix="#",
            has_header=False,
        )
        # If main chromosomes, plot.
        df_all = df_x.join(df_y, on="chrom", how="inner").filter(
            pl.col("chrom").str.contains_any(["MATERNAL", "PATERNAL"])
        )
        if args.omit_acrocentrics:
            df_all = df_all.filter(
                ~pl.col("chrom").str.contains_any(
                    ["chr13", "chr14", "chr15", "chr21", "chr22"]
                )
            )

        p = sns.jointplot(
            df_all,
            # df_other_chrs,
            x="qv_x",
            y="qv_y",
            kind="reg",
            marker="o",
            color=color,
            line_kws=dict(color=color, linestyle="dotted"),
        )
        # Highlight how acrocentric chromosomes are the outliers
        if not args.omit_acrocentrics:
            df_acro_chrs = df_all.filter(
                pl.col("chrom").str.contains_any(
                    ["chr13", "chr14", "chr15", "chr21", "chr22"]
                )
            )

            p.ax_joint.scatter(
                df_acro_chrs["qv_x"],
                df_acro_chrs["qv_y"],
                marker="o",
                facecolors="none",
                edgecolors="black",
            )

        ax = plt.gca()
        annotate(ax, df_all)
        ax.set_title(label, color=color, x=0.1, y=0.95, weight="bold")
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel)
        # annotate(ax, df_other_chrs)
        figs.append(p)

    fig = plt.figure(figsize=(20, 8))
    gs = gridspec.GridSpec(1, len(args.labels))

    for i, f in enumerate(figs):
        _ = SeabornFig2Grid(f, fig, gs[i])

    if not args.omit_acrocentrics:
        fig.legend(
            labels=["Acrocentric"],
            handles=[
                Line2D(
                    [],
                    [],
                    color="white",
                    marker="o",
                    markerfacecolor="none",
                    markeredgecolor="black",
                )
            ],
        )
    gs.tight_layout(fig)
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
