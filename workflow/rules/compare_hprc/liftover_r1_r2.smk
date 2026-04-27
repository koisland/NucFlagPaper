
rule align_r1_to_r2:
    input:
        ref=join(OUTPUT_DIR, "data", "{sm}_R2.fa.gz"),
        query=join(OUTPUT_DIR, "data", "{sm}_R1.fa.gz"),
    output:
        bam=join(OUTPUT_DIR, "data", "{sm}_R1_to_R2.bam"),
        paf=join(OUTPUT_DIR, "data", "{sm}_R1_to_R2.paf"),
    threads: 8
    conda:
        "../asm-to-reference-alignment/workflow/envs/env.yml"
    shell:
        """
        minimap2 -t {threads} -a --eqx --cs \
            -x asm5 --secondary=no -K 8G \
            {input.ref} {input.query} | \
        samtools view -F 4 -b - > {output.bam}
        samtools view -h {output.bam} | paftools.js sam2paf - > {output.paf}
        """


rule create_r1_to_r2_chain:
    input:
        rules.align_r1_to_r2.output.paf,
    output:
        join(OUTPUT_DIR, "data", "{sm}_R1_to_R2.chain"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        paf2chain -i {input} > {output}
        """


rule liftover_annotations_r1_to_r2:
    input:
        chain=rules.create_r1_to_r2_chain.output,
        segdups=expand(rules.download_segdup_annot.output, sm="{sm}", release="R2"),
    output:
        segdups=join(OUTPUT_DIR, "data", "{sm}_R1_segdups.bed"),
        segdups_unmapped=join(OUTPUT_DIR, "data", "{sm}_R1_segdups_unmapped.bed"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        liftOver {input.segdups} {input.chain} {output.segdups} {output.segdups_unmapped} -bedPlus=9
        """
