ASM_COLORS = {
    "CHM13_hifiasm_0.25.0": "red",
    "CHM13_verkko_2.2.1": "gray",
    "CHM13_verkko_2.3": "black",
}


rule draw_ng50:
    input:
        asm=[
            (
                expand(rules.denovo_verkko_output.output, asm=sm.split("_")[1], sm=sm)[
                    0
                ]
                + ".fai"
                if sm.split("_")[1] == "verkko"
                else expand(
                    rules.denovo_hifiasm_output.output, asm=sm.split("_")[1], sm=sm
                )[0]
                + ".fai"
            )
            for sm in config["samples"].keys()
        ],
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
        python {params.script} -i {input.asm} -l {params.labels} -c {params.colors} -o {output.plot}
        """
