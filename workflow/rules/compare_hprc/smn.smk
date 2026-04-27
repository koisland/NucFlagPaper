from typing import Iterator


# SMN2 and SMN1 (1 Mbp slop)
df_r1_segdups = pl.read_csv(config["annotations_segdups"]["R1"], separator="\t").select(
    sample=pl.col("sample"),
    haplotype=pl.col("haplotype"),
    location=pl.col("file_location"),
    release=pl.lit("R1"),
)
df_r2_segdups = pl.read_csv(config["annotations_segdups"]["R2"]).select(
    sample=pl.col("sample_id"),
    haplotype=pl.when(pl.col("haplotype") == 1)
    .then(pl.lit("paternal"))
    .otherwise(pl.lit("maternal")),
    location=pl.col("location"),
    release=pl.lit("R2"),
)
df_all_segdups = pl.concat([df_r1_segdups, df_r2_segdups])

# From other part of workflow
ASM_CHM13 = join(OUTPUT_DIR, "..", "assembly", "reference", "chm13v2.0.fa.gz")
BED_CHM13_SMN = (join("data", "compare_hprc", "chm13v2.0_smn_1mbp_slop.bed"),)
SM_METADATA = (join("data", "compare_hprc", "hprc_metadata.tsv"),)


def get_sms_release(*releases: str) -> Iterator[str]:
    yield from (f"{sm}_{release}" for sm in DATA.keys() for release in releases)


cmp_hprc_aln_cfg = {
    "ref": {"CHM13v2.0": ASM_CHM13},
    "sm": {
        f"{sm_release}": join(OUTPUT_DIR, "data", f"{sm_release}.fa.gz")
        for sm_release in get_sms_release(*RELEASES)
    },
    "temp_dir": join(OUTPUT_DIR, "aln", "temp"),
    "output_dir": join(OUTPUT_DIR, "aln"),
    "logs_dir": join(LOG_DIR, "aln"),
    "benchmarks_dir": join(BENCHMARK_DIR, "aln"),
    "aln_threads": 8,
    "aln_mem": "50GB",
    "mm2_opts": "-x asm20 --secondary=no -K 8G",
}


module align_asm_to_ref_hprc_r1_2:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        cmp_hprc_aln_cfg


use rule * from align_asm_to_ref_hprc_r1_2 as asm_ref_hprc_r1_r2_*


rule download_segdup_annot:
    output:
        join(OUTPUT_DIR, "data", "{sm}_{release}_sedef.bed"),
    params:
        uris=lambda wc: df_all_segdups.filter(
            pl.col("sample").eq(pl.lit(wc.sm))
            & pl.col("release").eq(pl.lit(wc.release))
        )["location"].to_list(),
    conda:
        # TODO: Consolidate envs
        "../../envs/curated.yaml"
    shell:
        """
        for uri in {params.uris}; do
            aws s3 cp ${{uri}} - >> {output}
        done
        """


rule query_w_impg_sm_to_chm13_paf:
    input:
        paf=lambda wc: expand(
            rules.asm_ref_hprc_r1_r2_bam_to_paf.output, ref=wc.ref, sm=wc.sm_release
        ),
        bed=BED_CHM13_SMN,
    output:
        join(OUTPUT_DIR, "smn", "{ref}", "{sm_release}.bed"),
    conda:
        "../compare_assembly/Snakemake-asm-evaluation-vg/workflow/envs/tools.yaml"
    shell:
        """
        impg query -p {input.paf} -b {input.bed} | \
        sort -k 1,1 -k2,2n | \
        bedtools merge -d 1000000 -i - -c 4,5,6,7,8,9,10 -o first,min,max,first,first,first,first > {output}
        """


include: "liftover_r1_r2.smk"


ruleorder: liftover_annotations_r1_to_r2 > download_segdup_annot


def get_segdup_annot(wc):
    sm, release = wc.sm_release.split("_", 1)
    # Cannot use segdup annotations from HPRC R1 because incomplete/incorrect.
    # Attempt to liftover from R2
    if release == "R2":
        return expand(rules.download_segdup_annot.output, sm=sm, release=release)
    else:
        return expand(rules.liftover_annotations_r1_to_r2.output.segdups, sm=sm)


rule intersect_sedef_nucflag:
    input:
        # Liftover of bed to new assembly coordinates
        bed=rules.query_w_impg_sm_to_chm13_paf.output,
        # Sedef annotations
        sedef=get_segdup_annot,
        # NucFlag output
        nucflag=lambda wc: expand(
            rules.run_nucflag.output,
            sm=wc.sm_release.split("_", 1)[0],
            release=wc.sm_release.split("_", 1)[1],
        ),
    output:
        sedef=join(OUTPUT_DIR, "smn", "{ref}", "{sm_release}_sedef_intersection.bed"),
        nucflag=join(
            OUTPUT_DIR, "smn", "{ref}", "{sm_release}_nucflag_intersection.bed"
        ),
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.sedef} -wb > {output.sedef}
        bedtools intersect -a {input.bed} -b {input.nucflag} -wb | \
            bedtools intersect -a - -b <(bedtools groupby -i {output.sedef} -g 1 -c 2,3 -o min,max) > {output.nucflag}
        """


# Looking for errors by length
# grep -v correct results/hprc/smn/CHM13v2.0/*R2_nucflag_intersection.bed | cut -f 1-3,14 | awk '{ print $0, $3-$2}' | sort -k5,5n


rule plot_r1_r2_chm13_coords:
    input:
        r1_sedef=expand(
            rules.intersect_sedef_nucflag.output.sedef,
            ref="CHM13v2.0",
            sm_release=list(get_sms_release("R1")),
        ),
        r1_nucflag=expand(
            rules.intersect_sedef_nucflag.output.nucflag,
            ref="CHM13v2.0",
            sm_release=list(get_sms_release("R1")),
        ),
        r2_sedef=expand(
            rules.intersect_sedef_nucflag.output.sedef,
            ref="CHM13v2.0",
            sm_release=list(get_sms_release("R2")),
        ),
        r2_nucflag=expand(
            rules.intersect_sedef_nucflag.output.nucflag,
            ref="CHM13v2.0",
            sm_release=list(get_sms_release("R2")),
        ),
        r1_fai=expand(
            join(OUTPUT_DIR, "data", "{sm_release}.fa.gz.fai"),
            sm_release=get_sms_release("R1"),
        ),
        r2_fai=expand(
            join(OUTPUT_DIR, "data", "{sm_release}.fa.gz.fai"),
            sm_release=get_sms_release("R2"),
        ),
        sample_metadata=SM_METADATA,
    output:
        multiext(join(OUTPUT_DIR, "smn", "smn_locus_r1_r2"), ".png", ".pdf"),
    params:
        script="workflow/scripts/compare_hprc/cmp_smn_r1_r2.py",
        output_prefix=lambda wc, output: splitext(output[0])[0],
    shell:
        """
        python {params.script} \
        --r1_sedef {input.r1_sedef} \
        --r2_sedef {input.r2_sedef} \
        --r1_nucflag {input.r1_nucflag} \
        --r2_nucflag {input.r2_nucflag} \
        --r1_fai {input.r1_fai} \
        --r2_fai {input.r2_fai} \
        --sm_metadata {input.sample_metadata} \
        --output_prefix {params.output_prefix}
        """


rule smn_all:
    input:
        rules.asm_ref_hprc_r1_r2_all.input,
        expand(
            rules.intersect_sedef_nucflag.output,
            ref="CHM13v2.0",
            sm_release=list(get_sms_release(*RELEASES)),
        ),
        rules.plot_r1_r2_chm13_coords.output,
