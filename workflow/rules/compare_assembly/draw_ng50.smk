

rule draw_ng50:
    input:
        fai=[get_assembly(sm) + ".fai" for sm in ASM_COLORS.keys()],
    output:
        plot=multiext(join(OUTPUT_DIR, "ng50", "plot_ng50"), ".png", ".pdf"),
    params:
        script="workflow/scripts/compare_assembly/generate_ng50.py",
        fai=lambda wc, input: " ".join(
            f"<(grep -v chrY {file})" if "chm13v2.0" in file else file
            for file in input.fai
        ),
        labels=" ".join(ASM_COLORS.keys()),
        colors=" ".join(ASM_COLORS.values()),
        output_prefix=lambda wc, output: splitext(output[0])[0],
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        python {params.script} -i {params.fai} -l {params.labels} -c {params.colors} -o {params.output_prefix}
        """
