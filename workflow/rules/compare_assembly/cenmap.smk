
rule run_cenmap:
    input:
        asm=lambda wc: (
            expand(
                rules.denovo_verkko_output.output, asm=wc.sm.split("_")[1], sm=wc.sm
            )[0]
            if wc.sm.split("_")[1] == "verkko"
            else expand(
                rules.denovo_hifiasm_output.output, asm=wc.sm.split("_")[1], sm=wc.sm
            )[0]
        ),
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


rule generate_simple_cytobands:
    input:
        rules.run_cenmap.output,
    output:
        join(OUTPUT_DIR, "cenmap", "{sm}_cytobands.bed"),
    shell:
        """
        # For each centromere, split at midpt into two annotations.
        awk -v OFS="\\t" '{{

            print $1, $2, $3, "acen"
        }}' {input}
        """


rule cenmap_all:
    input:
        expand(rules.run_cenmap.output, sm=config["samples"].keys()),
