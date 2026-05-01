import sys
import random
import argparse
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    loc="center left",
    alignment="left",
)
ANNOT_COLS = {
    "column_9": "qstrand",
    "column_11": "qchrom",
    "column_12": "qst",
    "column_13": "qend",
    "column_14": "name",
    "column_16": "strand",
    "column_19": "item_rgb",
    "column_4": "rchrom",
    "column_5": "rst",
    "column_6": "rend",
}
DUPMASKER_COLORS = {
    "chr1": "#ffe4e1",
    "chr2": "#000000",
    "chr3": "#2f4f4f",
    "chr4": "#191970",
    "chr5": "#6495ed",
    "chr6": "#0000cd",
    "chr7": "#00bfff",
    "chr8": "#4682b4",
    "chr9": "#00ffff",
    "chr10": "#e0ffff",
    "chr11": "#013220",
    "chr12": "#00ff7f",
    "chr13": "#7cfc00",
    "chr14": "#ffff00",
    "chr15": "#ffd700",
    "chr16": "#bc8f8f",
    "chr17": "#cd5c5c",
    "chr18": "#a52a2a",
    "chr19": "#ff6347",
    "chr20": "#cd853f",
    "chr21": "#ff0000",
    "chr22": "#ff69b4",
    "chrX": "#ff1493",
    "chrY": "#800000",
    "chr_random": "#9400d3",
}
# https://humanpangenome.org/samples/
AFR_POP_LABELS = {"ACB", "GWD", "ESN", "MKK", "LWK", "ASL", "YRI", "MSL", "ASW"}


def item_rgb_to_hex(color: str) -> str:
    return to_hex(tuple(int(c) / 255.0 for c in color.split(",")))


plt.rcParams["font.family"] = "Arial"


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False) -> None:
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
            labelbottom=False,
        )


def make_reorient_relative_df(df: pl.DataFrame) -> pl.DataFrame:
    """
    Reorient intervals based on CHM13 alignment.
    Also make relative. Requires all annotation data to be present as takes global minimum start position.
    """
    df_out = (
        df.with_columns(
            qst=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qend"))
            .otherwise(pl.col("qst")),
            qend=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qst"))
            .otherwise(pl.col("qend")),
            strand=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(
                pl.when(pl.col("strand").eq(pl.lit("-")))
                .then(pl.lit("+"))
                .when(pl.col("strand").eq(pl.lit("+")))
                .then(pl.lit("-"))
                .otherwise(pl.col("strand"))
            )
            .otherwise(pl.col("strand")),
        )
        .with_columns(
            st=(pl.col("qst") - pl.col("qst").min().over("qchrom")) + pl.col("rst"),
            end=(pl.col("qend") - pl.col("qst").min().over("qchrom")) + pl.col("rst"),
            mtch=pl.col("qchrom").str.extract_groups(r"^(?<sm>.*?)#(?<hap>1|2)#"),
        )
        .unnest("mtch")
        .sort(["rchrom", "hap", "st"])
    )
    # # Debugging
    # df_out_sm = df_out.filter(pl.col("qchrom").str.contains("HG01243#2"))
    # fig, ax = plt.subplots()
    # ax: Axes

    # for row in df_out_sm.iter_rows(named=True):
    #     ax.axvspan(row["st"], row["end"], color=item_rgb_to_hex(row["item_rgb"]))
    # fig.savefig("test.png")
    # df_out_sm.write_csv("test.csv")
    # breakpoint()
    return df_out


def draw_r1_r2_smn(
    df_all: pl.DataFrame, output_prefix: str, figsize: tuple[int, int] = (12, 24)
):
    sm_haps = sorted(set(df_all.select("sm", "hap").iter_rows()))
    fig, axes = plt.subplots(
        nrows=len(sm_haps) * 2,
        ncols=2,
        layout="constrained",
        sharex=True,
        height_ratios=[0.25, 1.0] * len(sm_haps),
        figsize=figsize,
    )
    for row_offset, annot in enumerate(("nucflag", "dupmasker")):
        df_annot = df_all.filter(pl.col("dtype").eq(pl.lit(annot)))
        for row, (sm, hap) in enumerate(sm_haps):
            row = (row * 2) + row_offset

            for col, release in enumerate(("R1", "R2")):
                df_sm_hap_annot = df_annot.filter(
                    pl.col("sm").eq(pl.lit(sm))
                    & pl.col("hap").eq(pl.lit(hap))
                    & pl.col("release").eq(pl.lit(release))
                )
                ax: Axes = axes[row, col]
                ax.set_ylim(0, 1)
                print(row, col, file=sys.stderr)

                ax.margins(x=0, y=0)

                minimalize_ax(ax, remove_ticks=True)

                polygons = []
                ht = ax.get_ylim()[1]
                for itv in df_sm_hap_annot.iter_rows(named=True):
                    hexcolor = item_rgb_to_hex(itv["item_rgb"])
                    if annot == "dupmasker":
                        # Draw triangle
                        # https://matplotlib.org/devdocs/api/_as_gen/matplotlib.patches.Polygon.html
                        if itv["strand"] == "+":
                            vertices = [
                                (itv["st"], 0),
                                (itv["st"], ht),
                                (itv["end"], ht / 2),
                            ]
                        else:
                            vertices = [
                                (itv["end"], 0),
                                (itv["end"], ht),
                                (itv["st"], ht / 2),
                            ]
                        poly = Polygon(xy=vertices, color=hexcolor)
                        polygons.append(poly)
                    else:
                        ax.axvspan(itv["st"], itv["end"], color=hexcolor)

                if annot == "dupmasker":
                    ax.add_collection(PatchCollection(polygons, match_original=True))

                if row_offset == 1:
                    has_break = False
                    for _, df_grp in df_sm_hap_annot.group_by(["qchrom"]):
                        # Is a break since if close to end of contig. Draw dashed line.
                        st_len_diff = abs(
                            df_grp["qst"].min() - df_grp["length"].first()
                        )
                        end_len_diff = abs(
                            df_grp["qend"].max() - df_grp["length"].first()
                        )
                        if st_len_diff < 50000:
                            x_line = df_grp["st"].min()
                        elif end_len_diff < 50000:
                            x_line = df_grp["end"].max()
                        else:
                            continue
                        # print(row, col, "break")
                        ax.axvline(
                            x=x_line,
                            ymin=0,
                            ymax=ht,
                            linestyle="dashed",
                            linewidth=1,
                            color="black",
                            zorder=2,
                        )
                        has_break = True
                    color = "red" if has_break else "black"
                    ax.set_ylabel(
                        f"{sm}_hap{hap}",
                        color=color,
                        rotation=0,
                        ha="right",
                        va="center",
                        fontsize=12,
                    )

    axes[0, 0].set_title("Release 1", fontsize=18)
    axes[0, 1].set_title("Release 2", fontsize=18)

    fig.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1_dupmasker", nargs="+", type=str)
    ap.add_argument("--r2_dupmasker", nargs="+", type=str)
    ap.add_argument("--r1_nucflag", nargs="+", type=str)
    ap.add_argument("--r2_nucflag", nargs="+", type=str)
    ap.add_argument("--r1_fai", nargs="+", type=str)
    ap.add_argument("--r2_fai", nargs="+", type=str)
    ap.add_argument("--n_subset", default=10, type=int)
    ap.add_argument("--seed", default=None, type=int)
    ap.add_argument("--output_prefix", default="out", type=str)
    args = ap.parse_args()

    dfs_r1_nucflag = []
    dfs_r2_nucflag = []
    dfs_r1_dupmasker = []
    dfs_r2_dupmasker = []

    # contig lengths
    all_lengths = []
    for lbl, fais in (("R1", args.r1_fai), ("R2", args.r2_fai)):
        for file in fais:
            df_fai = pl.read_csv(
                file,
                separator="\t",
                columns=[0, 1],
                new_columns=["chrom", "length"],
                has_header=False,
            )
            for chrom, length in df_fai.iter_rows():
                all_lengths.append([lbl, chrom, length])

    df_all_lengths = pl.DataFrame(
        all_lengths, schema=["release", "qchrom", "length"], orient="row"
    )

    for file in args.r1_nucflag:
        df_r1_nucflag = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
            )
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r1_nucflag.append(df_r1_nucflag)

    for file in args.r2_nucflag:
        df_r2_nucflag = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
            )
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r2_nucflag.append(df_r2_nucflag)

    for file in args.r1_dupmasker:
        df_r1_dupmasker = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r1_dupmasker.append(df_r1_dupmasker)

    for file in args.r2_dupmasker:
        df_r2_dupmasker = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(ANNOT_COLS.keys())
            .rename(ANNOT_COLS)
        )
        dfs_r2_dupmasker.append(df_r2_dupmasker)

    # Add length
    df_r1_nucflag = (
        pl.concat(dfs_r1_nucflag)
        .filter(pl.col("name").ne(pl.lit("correct")))
        .join(
            df_all_lengths.filter(pl.col("release").eq(pl.lit("R1"))),
            on="qchrom",
            how="left",
        )
    )
    df_r2_nucflag = (
        pl.concat(dfs_r2_nucflag)
        .filter(pl.col("name").ne(pl.lit("correct")))
        .join(
            df_all_lengths.filter(pl.col("release").eq(pl.lit("R1"))),
            on="qchrom",
            how="left",
        )
    )
    df_r1_dupmasker = pl.concat(dfs_r1_dupmasker).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("R1"))),
        on="qchrom",
        how="left",
    )
    df_r2_dupmasker = pl.concat(dfs_r2_dupmasker).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("R2"))),
        on="qchrom",
        how="left",
    )

    # Take min and max of same interval in case of alignment gap.
    # Also reorient if not forward oriented
    df_all = pl.concat(
        [
            df_r1_nucflag.with_columns(release=pl.lit("R1"), dtype=pl.lit("nucflag")),
            df_r2_nucflag.with_columns(release=pl.lit("R2"), dtype=pl.lit("nucflag")),
            df_r1_dupmasker.with_columns(
                release=pl.lit("R1"), dtype=pl.lit("dupmasker")
            ),
            df_r2_dupmasker.with_columns(
                release=pl.lit("R2"), dtype=pl.lit("dupmasker")
            ),
        ]
    )
    df_all = make_reorient_relative_df(df_all)

    nucflag_colors = {
        name: item_rgb_to_hex(color)
        for name, color in df_all.filter(pl.col("dtype").eq("nucflag"))
        .select("name", "item_rgb")
        .unique()
        .iter_rows()
    }

    # Draw subset of 10
    random.seed(args.seed)
    subset_sm = random.sample(list(df_all["sm"].unique()), args.n_subset)
    df_subset_all = df_all.filter(pl.col("sm").is_in(subset_sm))

    # draw_r1_r2_smn(df_all=df_all, output_prefix=args.output_prefix)
    draw_r1_r2_smn(
        df_all=df_subset_all,
        output_prefix=f"{args.output_prefix}_subset{args.n_subset}",
        figsize=(12, 8),
    )

    fig_legend, axes_legend = plt.subplots(nrows=3, ncols=1, layout="tight")
    axes_legend: list[Axes]
    for ax in axes_legend:
        minimalize_ax(ax, remove_ticks=True)
    axes_legend[0].legend(
        handles=[
            Patch(facecolor=color, edgecolor="black")
            for lbl, color in nucflag_colors.items()
        ],
        labels=nucflag_colors.keys(),
        ncols=len(nucflag_colors) // 2,
        title="NucFlag calls",
        **LEGEND_KWARGS,
    )
    axes_legend[1].legend(
        handles=[
            Patch(facecolor=color, edgecolor="black")
            for lbl, color in DUPMASKER_COLORS.items()
        ],
        labels=DUPMASKER_COLORS.keys(),
        ncols=len(DUPMASKER_COLORS) // 2,
        title="Segmental duplications",
        **LEGEND_KWARGS,
    )
    # https://www.geeksforgeeks.org/python/custom-legends-with-matplotlib/
    axes_legend[2].legend(
        handles=[Line2D([0], [0], color="black", linestyle="dashed")],
        labels=["Contig break"],
        **LEGEND_KWARGS | {"handlelength": 2.0},
    )
    fig_legend.savefig(f"{args.output_prefix}_legend.png", bbox_inches="tight", dpi=300)
    fig_legend.savefig(f"{args.output_prefix}_legend.pdf", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
