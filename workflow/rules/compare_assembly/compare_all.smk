

rule plot_compare_misassemblies:
    input:
        calls=[
            expand(
                rules.nf_denovo_check_asm_nucflag.output.misassemblies,
                sm=f"{sm}_hifi",
            )
            for sm in config["samples"].keys()
        ],
    output:
        plot=join(OUTPUT_DIR, "nucflag", "all_calls_compared.png"),
    params:
        script="workflow/scripts/compare_assembly/cmp_all_calls.py",
        labels=" ".join(config["samples"].keys()),
        colors=" ".join(ASM_COLORS.get(sm) for sm in config["samples"].keys()),
    shell:
        """
        python {params.script} -i {input.calls} -l {params.labels} -c {params.colors} -o {output.plot}
        """


rule overlap_annotation_data:
    input:
        censat=rules.liftover_annotations.output.censat,
        calls=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm="{sm}_hifi",
        ),
    output:
        censat=join(OUTPUT_DIR, "sequence_ctx", "{sm}_censat.bed"),
    shell:
        """
        bedtools intersect -loj -a <(cut -f1-4 {input.calls}) -b <(awk -v OFS="\\t" '{{ print $1, $2, $3, $4, $9 }}' {input.censat}) | \
            awk -v OFS="\\t" '{{ if ($5 == ".") {{ $8="no_overlap"; $9="128,128,128"; }}; print}}' > {output.censat}
        """


# TODO: Also generate context?


rule compare_all:
    input:
        rules.plot_compare_misassemblies.output,
