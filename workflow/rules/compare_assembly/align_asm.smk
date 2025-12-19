ALN_CFG = {
    "ref": {
        "CHM13_verkko": expand(
            rules.denovo_verkko_output.output, asm="verkko", sm="CHM13_verkko"
        )
    },
    "sm": {
        "CHM13_hifiasm": expand(
            rules.denovo_hifiasm_output.output, asm="hifiasm", sm="CHM13_hifiasm"
        )
    },
    "temp_dir": join(OUTPUT_DIR, "temp"),
    "output_dir": join(OUTPUT_DIR, "asm_ref"),
    "logs_dir": join(LOG_DIR, "asm_ref"),
    "benchmarks_dir": join(BENCHMARK_DIR, "asm_ref"),
    "aln_threads": 12,
    "aln_mem": "100GB",
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}


# Align hifiasm to verkko.
module align_asm_to_ref:
    snakefile:
        "asm-to-reference-alignment/workflow/Snakefile"
    config:
        ALN_CFG


use rule * from align_asm_to_ref as asm_ref_*


# https://dwinter.github.io/pafr/index.html
rule plot_regions:
    input:
        paf=[],
        regions=[],
        script="workflow/scripts/compare_assembly/plot_region_aln.R",
    output:
        dplot=[],
    params:
        bp_slop=1000,
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        Rscript {input.script} {input.paf} {input.regions} {output.dplot} {params.slop}
        """


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        rules.asm_ref_all.input,
    default_target: True
