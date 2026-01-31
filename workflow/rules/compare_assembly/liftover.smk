
wildcard_constraints:
    ref="CHM13v2.0",


rule download_annotations:
    output:
        gff=join(OUTPUT_DIR, "annot", "chm13.draft_v2.0.gene_annotation.gff3"),
        cytobands=join(OUTPUT_DIR, "annot", "chm13v2.0_cytobands_allchrs.bed"),
        censat=join(OUTPUT_DIR, "annot", "chm13v2.0_censat_v2.1.bed"),
        segdups=join(OUTPUT_DIR, "annot", "segDups_2024.bed"),
        segdups_bb=join(OUTPUT_DIR, "annot", "segDups_2024.bb"),
    params:
        url_gff="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3",
        url_cytobands="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed",
        url_censat="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed",
        url_sedef="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/segDups_2024.bb",
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        wget {params.url_gff} -O {output.gff}
        wget {params.url_cytobands} -O {output.cytobands}
        wget {params.url_censat} -O {output.censat}
        wget {params.url_sedef} -O {output.segdups_bb}
        bigbedtobed {output.segdups_bb} {output.segdups}
        """


saffire_cfg = {
    "ref": {"CHM13v2.0": rules.download_chm13.output.fa},
    "sm": {sm: get_assembly(sm) for sm in config["samples"].keys()},
    "temp_dir": join(OUTPUT_DIR, "saffire", "temp"),
    "output_dir": join(OUTPUT_DIR, "saffire"),
    "logs_dir": join(LOG_DIR, "saffire"),
    "benchmarks_dir": join(BENCHMARK_DIR, "saffire"),
    "aln_threads": 8,
    "aln_mem": "20GB",
    "mm2_opts": "-x asm5 --secondary=no -K 8G",
}


module align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        saffire_cfg


use rule * from align_asm_to_ref as asm_ref_*


rule paf2chain:
    input:
        rules.asm_ref_bam_to_paf.output,
    output:
        join(OUTPUT_DIR, "saffire", "{ref}", "chain", "{sm}.chain"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        paf2chain -i {input} > {output}
        """


rule liftover_annotations:
    input:
        chain=join(OUTPUT_DIR, "saffire", "{ref}", "chain", "{sm}.chain"),
        gff=rules.download_annotations.output.gff,
        cytobands=rules.download_annotations.output.cytobands,
        censat=rules.download_annotations.output.censat,
        segdups=rules.download_annotations.output.segdups,
    output:
        gff=join(OUTPUT_DIR, "annot", "{ref}_{sm}.gff"),
        gff_unmapped=join(OUTPUT_DIR, "annot", "{ref}_{sm}_unmapped.gff"),
        cytobands=join(OUTPUT_DIR, "annot", "{ref}_{sm}_cytobands.bed"),
        cytobands_unmapped=join(
            OUTPUT_DIR, "annot", "{ref}_{sm}_cytobands_unmapped.bed"
        ),
        censat=join(OUTPUT_DIR, "annot", "{ref}_{sm}_censat.bed"),
        censat_unmapped=join(OUTPUT_DIR, "annot", "{ref}_{sm}_censat_unmapped.bed"),
        segdups=join(OUTPUT_DIR, "annot", "{ref}_{sm}_segdups.bed"),
        segdups_unmapped=join(OUTPUT_DIR, "annot", "{ref}_{sm}_segdups_unmapped.bed"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        liftOver {input.gff} {input.chain} {output.gff} {output.gff_unmapped} -gff
        liftOver {input.cytobands} {input.chain} {output.cytobands} {output.cytobands_unmapped}
        liftOver {input.censat} {input.chain} {output.censat} {output.censat_unmapped}
        liftOver {input.segdups} {input.chain} {output.segdups} {output.segdups_unmapped} -bedPlus=9
        """


rule liftover_all:
    input:
        expand(
            rules.asm_ref_SafFire.output.bed,
            ref="CHM13v2.0",
            sm=config["samples"].keys(),
        ),
        rules.download_annotations.output,
        expand(rules.paf2chain.output, ref="CHM13v2.0", sm=config["samples"].keys()),
        expand(
            rules.liftover_annotations.output,
            ref="CHM13v2.0",
            sm=config["samples"].keys(),
        ),
