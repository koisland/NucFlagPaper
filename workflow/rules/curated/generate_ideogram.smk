# Download data source
rule download_annotation_data:
    output:
        segdups=join(IDEOGRAM_OUTPUT_DIR, "HG002.SDs.010624.45col.bb"),
        censat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.0.1.cenSatv2.0.noheader.bb"),
        cytoband_pat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.1.pat_cytoBandMapped.bb"),
        cytoband_mat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.1.mat_cytoBandMapped.bb"),
    params:
        output_dir=IDEOGRAM_OUTPUT_DIR,
        urls=" ".join(ANNOTATIONS_V101),
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        for url in {params.urls}; do
            fname=$(basename "${{url}}")
            new_path="{params.output_dir}/${{fname}}"
            wget "${{url}}" -O "${{new_path}}"
        done
        """


rule convert_annotations_to_bed:
    input:
        segdups=rules.download_annotation_data.output.segdups,
        censat=rules.download_annotation_data.output.censat,
        cytoband_pat=rules.download_annotation_data.output.cytoband_pat,
        cytoband_mat=rules.download_annotation_data.output.cytoband_mat,
    output:
        segdups=join(IDEOGRAM_OUTPUT_DIR, "HG002.SDs.010624.45col.bed"),
        censat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.0.1.cenSatv2.0.noheader.bed"),
        cytoband=join(IDEOGRAM_OUTPUT_DIR, "cytoBand.hg002v1.1.bed"),
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        bigbedtobed {input.segdups} {output.segdups}
        bigbedtobed {input.censat} {output.censat}
        bigbedtobed {input.cytoband_pat} {output.cytoband}.1
        bigbedtobed {input.cytoband_mat} {output.cytoband}.2
        cat {output.cytoband}.1 {output.cytoband}.2 > {output.cytoband}
        rm -f {output.cytoband}.1 {output.cytoband}.2
        """


saffire_cfg = {
    "ref": {
        "CHM13v2.0": ASM_CHM13
    },
    "sm": {
        "HG002v1.0.1": expand(rules.download_curated_asm.output.fa, version="v1.0.1")
    },
    "temp_dir": join(IDEOGRAM_OUTPUT_DIR, "saffire", "temp"),
    "output_dir": join(IDEOGRAM_OUTPUT_DIR, "saffire"),
    "logs_dir": join(IDEOGRAM_OUTPUT_DIR, "logs", "saffire"),
    "benchmarks_dir": join(IDEOGRAM_OUTPUT_DIR, "benchmarks", "saffire"),
    "aln_threads": 8,
    "aln_mem": "60GB",
    "mm2_opts": "-x asm20 --secondary=no -K 8G",
}


module liftover_align_asm_to_ref:
    snakefile:
        "../asm-to-reference-alignment/workflow/Snakefile"
    config:
        saffire_cfg


use rule * from liftover_align_asm_to_ref as liftover_from_chm13_*


rule liftover_from_chm13_paf2chain:
    input:
        rules.liftover_from_chm13_bam_to_paf.output,
    output:
        join(IDEOGRAM_OUTPUT_DIR, "saffire", "{ref}", "chain", "{sm}.chain"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        paf2chain -i {input} > {output}
        """


# Cytobands for HG002 v1.0.1 are incomplete
# Cytobands for HG002 v1.1 chrY are likely incorrect in Yq12.
# We liftover and manually correct the Yq12
rule liftover_cytobands_from_chm13:
    input:
        chain=rules.liftover_from_chm13_paf2chain.output,
        cytobands_chm13=CYTOBANDS_CHM13,
    output:
        cytobands=join(IDEOGRAM_OUTPUT_DIR, "{ref}_{sm}_cytobands.bed"),
        cytobands_unmapped=join(
            IDEOGRAM_OUTPUT_DIR, "{ref}_{sm}_cytobands_unmapped.bed"
        ),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        liftOver -multiple {input.cytobands_chm13} {input.chain} {output.cytobands}.tmp {output.cytobands_unmapped}
        sort -k1,1 -k2,2n {output.cytobands}.tmp > {output.cytobands}
        rm -f {output.cytobands}.tmp
        """


rule plot_ideogram:
    input:
        script="workflow/scripts/metrics/plot_curated_ideogram.py",
        cytobands=rules.convert_annotations_to_bed.output.cytoband,
        curated=expand(rules.convert_vcf_to_bed.output, version="v1.0.1"),
        nucflag_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="nucflag",
            version="v1.0.1",
        ),
        inspector_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="inspector",
            version="v1.0.1",
        ),
        flagger_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="flagger",
            version="v1.0.1",
        ),
        deepvariant_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="deepvariant",
            version="v1.0.1",
        ),
        fai=expand(rules.download_curated_asm.output.fai, version=f"v1.0.1")[0],
    output:
        plots=[
            join(IDEOGRAM_OUTPUT_DIR, "ideogram_{cond}.png"),
            join(IDEOGRAM_OUTPUT_DIR, "ideogram_{cond}_overview.png"),
            join(IDEOGRAM_OUTPUT_DIR, "ideogram_{cond}.pdf"),
        ],
    conda:
        "../../envs/tools.yaml"
    params:
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
        include_false_positives=lambda wc: "-f" if wc.cond == "fp" else "",
    shell:
        """
        python {input.script} \
        --fai {input.fai} \
        --cytobands {input.cytobands} \
        --curated {input.curated} \
        --nucflag {input.nucflag_calls}/missed_calls.bed \
        --inspector {input.inspector_calls}/missed_calls.bed \
        --deepvariant {input.deepvariant_calls}/missed_calls.bed \
        --flagger {input.flagger_calls}/missed_calls.bed \
        {params.include_false_positives} \
        -o {params.output_prefix}
        """


rule plot_ideogram_chrom:
    input:
        script="workflow/scripts/metrics/plot_curated_ideogram_chrom.py",
        cytobands=rules.convert_annotations_to_bed.output.cytoband,
        curated=expand(rules.convert_vcf_to_bed.output, version=f"v1.0.1"),
        nucflag_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="nucflag",
            version="v1.0.1",
        ),
        inspector_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="inspector",
            version="v1.0.1",
        ),
        flagger_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="flagger",
            version="v1.0.1",
        ),
        deepvariant_calls=expand(
            rules.calculate_precision_recall.output.missed_calls_dir,
            tool="deepvariant",
            version="v1.0.1",
        ),
        segdups=rules.convert_annotations_to_bed.output.segdups,
        censat=rules.convert_annotations_to_bed.output.censat,
    wildcard_constraints:
        chrom="|".join(CHROMS),
    output:
        png=join(IDEOGRAM_OUTPUT_DIR, "{chrom}_{cond}.png"),
        pdf=join(IDEOGRAM_OUTPUT_DIR, "{chrom}_{cond}.pdf"),
    params:
        include_false_positives=lambda wc: "-f" if wc.cond == "fp" else "",
        output_prefix=lambda wc, output: splitext(output[0])[0],
        # Set to length of chrY
        xlim=lambda wc: (
            "--xlim '(0, 62425491)'"
            if wc.chrom == "chr21_PATERNAL" or wc.chrom == "chrY_PATERNAL"
            else ""
        ),
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        python {input.script} \
        --cytobands {input.cytobands} \
        --curated {input.curated} \
        --nucflag {input.nucflag_calls}/missed_calls.bed \
        --inspector {input.inspector_calls}/missed_calls.bed \
        --flagger {input.flagger_calls}/missed_calls.bed \
        --deepvariant {input.deepvariant_calls}/missed_calls.bed \
        --segdups {input.segdups} \
        --censat {input.censat} \
        -o {params.output_prefix} {params.include_false_positives} \
        -c {wildcards.chrom} {params.xlim}
        """


rule generate_ideogram:
    input:
        expand(rules.plot_ideogram.output, cond=["fp", "none"]),
        expand(
            rules.plot_ideogram_chrom.output,
            chrom=[chrom for chrom in CHROMS if chrom != "all"],
            cond=["fp", "none"],
        ),
        # Liftover CHM13 cytobands to fix chrY Yq12 manually.
        expand(
            rules.liftover_cytobands_from_chm13.output,
            ref="CHM13v2.0",
            sm="HG002v1.0.1",
        ),
    default_target: True
