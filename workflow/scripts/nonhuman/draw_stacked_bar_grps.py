import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch

LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    alignment="left",
    edgecolor="black",
)
LARGE_ERRORS = (
    "false_dup",
    "misjoin",
    "collapse",
    "scaffold",
    "softclip",
)
SMALL_ERRORS = (
    "insertion",
    "deletion",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "other_repeat",
    "mismatch",
    "low_quality",
    "het_or_mismap",
)
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["text.usetex"] = False


def draw_grouped_bar_w_calls_by_asm(
    df: pl.DataFrame,
    asm_palette: dict[str, str],
    call_palette: dict[str, str],
    output_prefix: str,
):
    df_grouped_calls = (
        df.unpivot(
            index=["asm", "group", "group_length", "group_rows"],
            variable_name="name",
            value_name="count",
        )
        .filter(pl.col("name").ne(pl.lit("correct")))
        .drop("group_length", "group_rows")
    )
    df_all_grouped_calls = df_grouped_calls.group_by(["asm", "name"]).agg(
        pl.col("count").sum()
    )
    fig_grouped_calls, ax_grouped_calls = plt.subplots(
        layout="constrained", figsize=(12, 3)
    )
    ax_grouped_calls: Axes
    # Draw cen bar
    bars_cens = sns.barplot(
        data=df_grouped_calls.filter(pl.col("group").eq(pl.lit("cen"))),
        x="name",
        y="count",
        hue="asm",
        palette=asm_palette,
        order=[*SMALL_ERRORS, *LARGE_ERRORS],
        hue_order=asm_palette.keys(),
        legend=None,
        ax=ax_grouped_calls,
        zorder=2,
    )
    for bar in bars_cens.patches:
        bar.set_hatch("////")

    sns.barplot(
        data=df_all_grouped_calls,
        x="name",
        y="count",
        hue="asm",
        order=[*SMALL_ERRORS, *LARGE_ERRORS],
        hue_order=asm_palette.keys(),
        palette=asm_palette,
        legend="full",
        ax=ax_grouped_calls,
    )
    for c in ax_grouped_calls.containers:
        labels = [f"{round(v.get_height())}" for v in c]
        ax_grouped_calls.bar_label(
            c,
            labels=labels,
            label_type="edge",
            fontsize=10,
        )

    ax_grouped_calls.set_xlabel(None)
    ax_grouped_calls.set_ylabel("# of calls", fontsize=14)
    for lbl in ax_grouped_calls.xaxis.get_majorticklabels():
        # Color and stroke around text
        lbl.set_color(call_palette[lbl.get_text()])
        lbl.set_fontsize(14)
        # lbl.set_path_effects([pe.Stroke(linewidth=0.1, foreground="black")])
        lbl.set_rotation(45)
        lbl.set_horizontalalignment("right")
        lbl.set_rotation_mode("anchor")

    for spine in ("top", "right"):
        ax_grouped_calls.spines[spine].set_visible(False)

    handles, labels = ax_grouped_calls.get_legend_handles_labels()
    for handle in handles:
        handle.set_edgecolor("black")
    handles.append(Patch(hatch="////", edgecolor="black", fill=False))
    labels.append("Centromere")
    ax_grouped_calls.legend(handles, labels)
    sns.move_legend(
        ax_grouped_calls, loc="upper right", fontsize=14, title=None, **LEGEND_KWARGS
    )
    fig_grouped_calls.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    fig_grouped_calls.savefig(f"{output_prefix}.pdf", dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", help="Input status groups summary.")
    ap.add_argument(
        "-k", "--color_key", help="Misassembly color key as two column TSV."
    )
    ap.add_argument("-l", "--labels", nargs="+", required=True)
    ap.add_argument("-c", "--colors", nargs="+", required=True)
    ap.add_argument("-o", "--output_prefix", help="Output prefix.")
    args = ap.parse_args()
    asm_colors = dict(zip(args.labels, args.colors))

    df_status = (
        pl.read_csv(
            args.infile,
            separator="\t",
        )
        .with_columns(
            asm=pl.col("file").str.replace(
                r"_(hifi|ont)_status_grps_(count|length)", ""
            )
        )
        .drop("file")
    )
    # Plot bar
    color_key = {
        name: rgb2hex([int(e) / 255.0 for e in itemRgb.split(",")])
        if not itemRgb.startswith("#")
        else itemRgb
        for name, itemRgb in pl.read_csv(
            args.color_key, has_header=False, separator="\t"
        ).iter_rows()
    }
    draw_grouped_bar_w_calls_by_asm(
        df_status,
        asm_palette=asm_colors,
        call_palette=color_key,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
