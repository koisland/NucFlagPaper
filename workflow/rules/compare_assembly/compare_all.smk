

rule plot_compare_misassemblies:
    input:
        calls=[
            expand(
                rules.nf_denovo_check_asm_nucflag.output.misassemblies,
                sm="{sm}_hifi",
            )
            for sm in config["samples"].keys()
        ],
    output:
        bed=join(OUTPUT_DIR, "nucflag", "all_calls_compared.png"),
    params:
        script="workflow/scripts/compare_assembly/cmp_all_calls.py",
        labels=" ".join(config["samples"].keys()),
        colors=" ".join(ASM_COLORS.get(sm) for sm in config["samples"].keys()),
    shell:
        """
        python {params.script} -i {input.calls} -l {params.labels} -c {params.colors}
        """
