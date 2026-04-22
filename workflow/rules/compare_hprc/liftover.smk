# SMN2 and SMN1 (1 Mbp slop)
# CHM13v2.0 chr5:69,359,653-72,931,727


# From other part of workflow
ASM_CHM13 = join(OUTPUT_DIR, "assembly", "reference", "chm13v2.0.fa.gz")

cmp_hprc_saffire_cfg = {
    "ref": {"CHM13v2.0": ASM_CHM13},
    "sm": {
        sm: join(OUTPUT_DIR, "data", f"{sm}_{release}.fa.gz")
        for sm in DATA.keys()
        for release in RELEASES
    },
    "temp_dir": join(OUTPUT_DIR, "saffire", "temp"),
    "output_dir": join(OUTPUT_DIR, "saffire"),
    "logs_dir": join(LOG_DIR, "saffire"),
    "benchmarks_dir": join(BENCHMARK_DIR, "saffire"),
    "aln_threads": 8,
    "aln_mem": "20GB",
    "mm2_opts": "-x asm20 --secondary=no -K 8G",
}


module align_asm_to_ref_hprc_r1_2:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        saffire_cfg


use rule * from align_asm_to_ref_hprc_r1_2 as asm_ref_hprc_r1_r2_*


rule query_w_impg_sm_to_chm13_paf:
    input:
        paf=rules.asm_ref_hprc_r1_r2_bam_to_paf.output,
    output:
        join(OUTPUT_DIR, "saffire", "{ref}", "chain", "{sm}.chain"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        impg query -p {input.paf}
        """
