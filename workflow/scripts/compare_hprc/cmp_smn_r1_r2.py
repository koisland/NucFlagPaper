import argparse
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.colors import to_hex


LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    loc="center left",
    alignment="left",
)
NUCFLAG_COLS = {
    "column_9": "qstrand",
    "column_11": "qchrom",
    "column_12": "qst",
    "column_13": "qend",
    "column_14": "name",
    "column_19": "item_rgb",
    "column_4": "rchrom",
    "column_5": "rst",
    "column_6": "rend",
}
SEDEF_COLS = {
    "column_9": "qstrand",
    "column_11": "qchrom",
    "column_12": "qst",
    "column_13": "qend",
    "column_34": "identity",
    "column_19": "item_rgb",
    "column_4": "rchrom",
    "column_5": "rst",
    "column_6": "rend",
}
# https://humanpangenome.org/samples/
AFR_POP_LABELS = {"ACB", "GWD", "ESN", "MKK", "LWK", "ASL", "YRI", "MSL", "ASW"}


def item_rgb_to_hex(color: str) -> str:
    return to_hex(tuple(int(c) / 255.0 for c in color.split(",")))


# From https://github.com/EichlerLab/sedef_smk/blob/46f009cbe319613cb0d1f9a26393e25bf1f2c309/scripts/sedef_to_bed.py#L31-L46
SEGDUP_LEGEND = {
    r"Greater than 99%": item_rgb_to_hex("255,103,0"),
    r"98% - 99%": item_rgb_to_hex("204,204,0"),
    # Midpoint between min and max identity
    r"90% - 98%": item_rgb_to_hex("162,162,162"),
    r"Less than 90%": item_rgb_to_hex("147,112,219"),
}

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


def reformat_nucflag_df(df: pl.DataFrame) -> pl.DataFrame:
    """
    Reorient intervals based on strand.
    Take minimum and maximum coordinates of reference
    """
    return (
        df.with_columns(
            qst=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qend"))
            .otherwise(pl.col("qst")),
            qend=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qst"))
            .otherwise(pl.col("qend")),
        )
        .group_by(["qchrom", "qst", "qend", "name", "item_rgb", "rchrom", "length"])
        .agg(pl.col("rst").min(), pl.col("rend").max())
        .with_columns(
            st=pl.col("qst") + pl.col("rst"),
            end=pl.col("qend") + pl.col("rst"),
            mtch=pl.col("qchrom").str.extract_groups(r"^(?<sm>.*?)#(?<hap>1|2)#"),
        )
        .unnest("mtch")
        .sort(["rchrom", "hap", "st"])
    )


def reformat_sedef_df(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.with_columns(
            qst=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qend"))
            .otherwise(pl.col("qst")),
            qend=pl.when(pl.col("qstrand").eq(pl.lit("-")))
            .then(pl.col("length") - pl.col("qst"))
            .otherwise(pl.col("qend")),
        )
        .group_by(["qchrom", "qst", "qend", "item_rgb", "rchrom", "length"])
        .agg(pl.col("rst").min(), pl.col("rend").max())
        .with_columns(
            st=pl.col("qst") + pl.col("rst"),
            end=pl.col("qend") + pl.col("rst"),
            mtch=pl.col("qchrom").str.extract_groups(r"^(?<sm>.*?)#(?<hap>1|2)#"),
        )
        .unnest("mtch")
        .sort(["rchrom", "hap", "st"])
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1_sedef", nargs="+", type=str)
    ap.add_argument("--r2_sedef", nargs="+", type=str)
    ap.add_argument("--r1_nucflag", nargs="+", type=str)
    ap.add_argument("--r2_nucflag", nargs="+", type=str)
    ap.add_argument("--r1_fai", nargs="+", type=str)
    ap.add_argument("--r2_fai", nargs="+", type=str)
    ap.add_argument("--sample_metadata", type=str)
    ap.add_argument("--output_prefix", default="out", type=str)
    args = ap.parse_args()

    dfs_r1_nucflag = []
    dfs_r2_nucflag = []
    dfs_r1_sedef = []
    dfs_r2_sedef = []

    # Population metadata
    samples_afr = set(
        pl.read_csv(args.sample_metadata, separator="\t")
        .filter(pl.col("Population Abbreviation").is_in(AFR_POP_LABELS))
        .get_column("Sample ID")
    )

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
            .select(NUCFLAG_COLS.keys())
            .rename(NUCFLAG_COLS)
        )
        dfs_r1_nucflag.append(df_r1_nucflag)

    for file in args.r2_nucflag:
        df_r2_nucflag = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
            )
            .select(NUCFLAG_COLS.keys())
            .rename(NUCFLAG_COLS)
        )
        dfs_r2_nucflag.append(df_r2_nucflag)

    for file in args.r1_sedef:
        df_r1_sedef = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(SEDEF_COLS.keys())
            .rename(SEDEF_COLS)
            .drop("identity")
        )
        dfs_r1_sedef.append(df_r1_sedef)

    for file in args.r2_sedef:
        df_r2_sedef = (
            pl.read_csv(file, separator="\t", has_header=False)
            .select(SEDEF_COLS.keys())
            .rename(SEDEF_COLS)
            .drop("identity")
        )
        dfs_r2_sedef.append(df_r2_sedef)

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
    df_r1_sedef = pl.concat(dfs_r1_sedef).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("R1"))),
        on="qchrom",
        how="left",
    )
    df_r2_sedef = pl.concat(dfs_r2_sedef).join(
        df_all_lengths.filter(pl.col("release").eq(pl.lit("R2"))),
        on="qchrom",
        how="left",
    )

    # Take min and max of same interval in case of alignment gap.
    # Also reorient if not forward oriented
    df_nucflag = pl.concat(
        [
            reformat_nucflag_df(df_r1_nucflag).with_columns(release=pl.lit("R1")),
            reformat_nucflag_df(df_r2_nucflag).with_columns(
                release=pl.lit("R2"),
            ),
        ]
    )
    df_sedef = pl.concat(
        [
            reformat_sedef_df(df_r1_sedef).with_columns(
                release=pl.lit("R1"),
            ),
            reformat_sedef_df(df_r2_sedef).with_columns(
                release=pl.lit("R2"),
            ),
        ]
    )

    df_min_st = (
        pl.concat((df_nucflag.select("qchrom", "st"), df_sedef.select("qchrom", "st")))
        .group_by(["qchrom"])
        .agg(min_st=pl.col("st").min())
    )
    # Convert to relative coordinates
    df_nucflag = df_nucflag.join(df_min_st, on="qchrom", how="left").with_columns(
        pl.col("st") - pl.col("min_st"), pl.col("end") - pl.col("min_st")
    )
    nucflag_colors = {
        name: item_rgb_to_hex(color)
        for name, color in df_nucflag.select("name", "item_rgb").unique().iter_rows()
    }
    df_sedef = df_sedef.join(df_min_st, on="qchrom", how="left").with_columns(
        pl.col("st") - pl.col("min_st"), pl.col("end") - pl.col("min_st")
    )

    sm_haps = sorted(set(df_nucflag.select("sm", "hap").iter_rows()))
    fig, axes = plt.subplots(
        nrows=len(sm_haps) * 2,
        ncols=2,
        layout="constrained",
        sharex=True,
        height_ratios=[0.25, 1.0] * len(sm_haps),
        figsize=(12, 36),
    )
    for row_offset, df_annot in enumerate((df_nucflag, df_sedef)):
        for row, (sm, hap) in enumerate(sm_haps):
            row = (row * 2) + row_offset
            # Color if sample of African ancestry
            label_color = "#E8A952" if sm in samples_afr else "black"

            for col, release in enumerate(("R1", "R2")):
                df_sm_hap_annot = df_annot.filter(
                    pl.col("sm").eq(pl.lit(sm))
                    & pl.col("hap").eq(pl.lit(hap))
                    & pl.col("release").eq(pl.lit(release))
                )
                ax: Axes = axes[row, col]
                # TODO: Add continental group
                if row_offset == 1:
                    ax.set_ylabel(
                        f"{sm}_hap{hap}",
                        rotation=0,
                        ha="right",
                        va="center",
                        color=label_color,
                        fontsize=12,
                    )

                ax.margins(x=0, y=0)

                minimalize_ax(ax, remove_ticks=True)

                for itv in df_sm_hap_annot.iter_rows(named=True):
                    hexcolor = item_rgb_to_hex(itv["item_rgb"])
                    ax.axvspan(itv["st"], itv["end"], color=hexcolor)

                if row_offset == 1:
                    for _, df_grp in df_sm_hap_annot.group_by(["qchrom"]):
                        # Is a break since if close to end of contig. Draw dashed line.
                        len_diff = abs(df_grp["qend"].max() - df_grp["length"].first())
                        if len_diff < 50000:
                            # print(row, col, "break")
                            max_end = df_grp["end"].max()
                            ax.axvline(
                                x=max_end,
                                linestyle="dashed",
                                linewidth=1,
                                color="black",
                                zorder=2,
                            )

    axes[0, 0].set_title("Release 1", fontsize=18)
    axes[0, 1].set_title("Release 2", fontsize=18)

    fig_legend, axes_legend = plt.subplots(nrows=2, ncols=1, layout="tight")
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
            for lbl, color in SEGDUP_LEGEND.items()
        ],
        labels=SEGDUP_LEGEND.keys(),
        ncols=len(SEGDUP_LEGEND),
        title="Segemental duplications",
        **LEGEND_KWARGS,
    )

    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=300)
    fig_legend.savefig(f"{args.output_prefix}_legend.png", bbox_inches="tight", dpi=300)
    fig_legend.savefig(f"{args.output_prefix}_legend.pdf", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    raise SystemExit(main())
