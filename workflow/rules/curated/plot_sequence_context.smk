
rule overlap_annotation_data:
    input:
        censat=rules.convert_to_bed.output.censat,
        segdup=rules.convert_to_bed.output.segdups,
        calls=get_bed_files,
    output:
        censat=join(OUTPUT_DIR, "sequence_ctx", "{tool}_{version}_censat.bed"),
        segdup=join(OUTPUT_DIR, "sequence_ctx", "{tool}_{version}_segdup.bed"),
    shell:
        """
        bedtools intersect -loj -a <(cut -f1-4 {input.calls}) -b <(awk -v OFS="\\t" '{{ print $1, $2, $3, $4, $9 }}' {input.censat}) | \
            awk -v OFS="\\t" '{{ if ($5 == ".") {{ $8="no_overlap"; $9="128,128,128"; }}; print}}' > {output.censat}
        bedtools intersect -loj -a <(cut -f1-4 {input.calls}) -b <(awk -v OFS="\\t" '{{ print $1, $2, $3, $37, $9 }}' {input.segdup}) | \
            awk -v OFS="\\t" '{{ if ($5 == ".") {{ $8="no_overlap"; $9="128,128,128"; }}; print}}' > {output.segdup}
        """


rule plot_seq_context:
    input:
        calls=get_bed_files,
        censat=rules.overlap_annotation_data.output.censat,
        segdup=rules.overlap_annotation_data.output.segdup,
    output:
        all_plot=join(OUTPUT_DIR, "sequence_ctx", "{tool}_{version}_total_{chrom}.png"),
        censat_plot=join(
            OUTPUT_DIR, "sequence_ctx", "{tool}_{version}_censat_{chrom}.png"
        ),
        segdup_plot=join(
            OUTPUT_DIR, "sequence_ctx", "{tool}_{version}_segdup_{chrom}.png"
        ),
    params:
        script="workflow/scripts/curated/plot_seq_context.py",
        chrom="{chrom}",
        output_prefix=lambda wc, output: join(
            dirname(output[0]), f"{wc.tool}_{wc.version}"
        ),
    shell:
        """
        python {params.script} \
        {params.chrom} \
        {input.calls} {input.censat} {input.segdup} \
        {params.output_prefix}
        """


rule sequence_context_all:
    input:
        expand(
            rules.plot_seq_context.output,
            version=["v1.0.1"],
            tool=["flagger", "inspector", "nucflag", "deepvariant"],
            chrom=CHROMS,
        ),
