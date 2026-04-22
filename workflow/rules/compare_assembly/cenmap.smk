
rule run_cenmap:
    input:
        asm=lambda wc: get_assembly(wc.sm),
    output:
        join(
            OUTPUT_DIR,
            "cenmap",
            "{sm}",
            "{sm}",
            "results",
            "final",
            "bed",
            "{sm}_complete_correct_cens.bed",
        ),
    params:
        output_dir=join(OUTPUT_DIR, "cenmap", "{sm}"),
        profile="workflow/profiles/lpc_all",
    conda:
        "../../envs/cenmap.yaml"
    log:
        join(LOG_DIR, "cenmap", "{sm}.log"),
    threads: 1
    shell:
        """
        cenmap \
        -i {input.asm} \
        -o {params.output_dir} \
        -j 20 \
        --workflow-profile {params.profile} \
        -s {wildcards.sm} 2> {log}
        """


rule cenmap_all:
    input:
        expand(rules.run_cenmap.output, sm=config["samples"].keys()),
