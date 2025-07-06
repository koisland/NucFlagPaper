import os
import sys
import glob
import argparse
import polars as pl

"""
python misasim/calculate_precision_recall.py \
-p /project/logsdon_shared/projects/NucFlagPaper/results/misasim/simulated/params.tsv \
-m /project/logsdon_shared/projects/NucFlagPaper/results/misasim/simulated \
-n /project/logsdon_shared/projects/NucFlagPaper/results/nucflag_misasim
"""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--params_file", required=True)
    ap.add_argument("-m", "--simulated_bed_dir", required=True)
    ap.add_argument("-n", "--nucflag_bed_dir", required=True)

    args = ap.parse_args()

    df_params = pl.read_csv(
        args.params_file,
        separator="\t"
    )

    for file in glob.glob(os.path.join(args.simulated_bed_dir, "*.bed")):
        seed = int(os.path.splitext(os.path.basename(file))[0])
        seed, mtype, num, max_length = df_params.filter(pl.col("seed") == seed).row(0)
        columns = (
            ("chrom", "st", "end", "num_dupes", "seq")
            if mtype == "false-duplication"
            else
            ("chrom", "st", "end", "seq")
        )
        df_misasim = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
                new_columns=columns
            )
            .with_columns(seed=seed)
            .cast({"seed": pl.Int64})
            .join(
                df_params,
                on=["seed"],
                how="left",
            )
        )
        chroms_misasim = df_misasim["chrom"].unique()

        nucflag_beds = glob.glob(
            os.path.join(args.nucflag_bed_dir, f"*-{mtype}-{num}-{max_length}_misassemblies.bed")
        )
        for bed in nucflag_beds:
            chrom_perc_mtype, bed_num, bed_max_length = os.path.basename(bed).replace("_misassemblies.bed", "").rsplit("-", 2)
            _, bed_perc, bed_mtype = chrom_perc_mtype.split("-", 2)

            df_nucflag = pl.read_csv(
                bed,
                separator="\t",
                has_header=False,
                new_columns=("chrom", "st", "end", "mtype")
            )
            for chrom, df_chrom_misasim in df_misasim.group_by(["chrom"]):
                chrom = chrom[0]
                print(bed_num, bed_max_length, bed_perc, bed_mtype)
                df_chrom_nucflag = df_nucflag.filter(pl.col("chrom") == chrom)

                breakpoint()

    # Find all files in -n with "_misassemblies.bed"

if __name__ == "__main__":
    raise SystemExit(main())
