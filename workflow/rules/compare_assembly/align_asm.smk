

rule merge_asm:
    input:
        asm_1="",
        asm_2="",
    output:
        merged_asms="",
    shell:
        """
        cat {input.asm_1} {input.asm_2} > {output}
        """


rule align_compare_asm:
    input:
        asm_1="",
        asm_2="",
    output:
        "",
    shell:
        """
        """
