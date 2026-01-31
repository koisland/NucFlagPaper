

rule draw_ng50:
    input:
        asm=[get_assembly(sm) + ".fai" for sm in config["samples"].keys()],
        ref=rules.download_chm13.output.fai,
    output:
        plot=join(OUTPUT_DIR, "ng50", "plot_ng50.png"),
    params:
        script="workflow/scripts/compare_assembly/generate_ng50.py",
        labels=" ".join(config["samples"].keys()),
        colors=" ".join(ASM_COLORS.get(sm, "None") for sm in config["samples"].keys()),
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        python {params.script} -i {input.asm} <(grep -v chrY {input.ref}) -l {params.labels} CHM13v2.0 -c {params.colors} red -o {output.plot}
        """
