import sys
import polars as pl

r1_all_manifest = "/project/logsdon_shared/projects/Keith/NucFlagPaper/config/hprc_r1_hifi_alignments.csv"
r2_asm_manifest = (
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/config/hprc_r2_assemblies.csv"
)
r2_hifi_mapping_manifest = "/project/logsdon_shared/projects/Keith/NucFlagPaper/config/hprc_r2_hifi_alignments.csv"

# https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables#outstanding-samples
df_all_r1 = (
    pl.read_csv(r1_all_manifest)
    .with_columns(release=pl.lit("R1"))
    .rename({"hifi_mapping": "bam"})
    .filter(pl.col("sample") != "HG03492")
)

df_asm_r2 = (
    pl.read_csv(r2_asm_manifest)
    .select("sample_id", "haplotype", "assembly")
    .pivot(on="haplotype", index="sample_id", values="assembly")
    .with_columns(
        pl.when(pl.col("1").is_null())
        .then(pl.col("0"))
        .otherwise(pl.col("1"))
        .alias("1")
    )
    .select("sample_id", "1", "2")
)
df_hifi_r2 = pl.read_csv(r2_hifi_mapping_manifest).select("sample_id", "bam")

df_all_r2 = (
    df_asm_r2.join(df_hifi_r2, how="left", on="sample_id")
    .with_columns(release=pl.lit("R2"))
    .rename({"1": "hap1_fasta", "2": "hap2_fasta", "sample_id": "sample"})
    .filter(pl.col("sample").is_in(df_all_r1["sample"]))
)

df_all = pl.concat([df_all_r1, df_all_r2])
df_all.write_csv(sys.stdout)
