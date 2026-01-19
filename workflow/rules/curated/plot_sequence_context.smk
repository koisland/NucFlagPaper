from snakemake.io import Wildcards


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
        calls=[
            get_bed_files(Wildcards(fromdict=dict(tool=tool, version="v1.0.1")))
            for tool in TOOLS.keys()
        ],
        censat=expand(
            rules.overlap_annotation_data.output.censat,
            tool=TOOLS.keys(),
            version="v1.0.1",
        ),
        segdup=expand(
            rules.overlap_annotation_data.output.segdup,
            tool=TOOLS.keys(),
            version="v1.0.1",
        ),
    output:
        all_plot=join(OUTPUT_DIR, "sequence_ctx", "{chrom}.png"),
        ctx_plot=join(OUTPUT_DIR, "sequence_ctx", "{chrom}_ctx.png"),
    params:
        script="workflow/scripts/curated/plot_seq_context.py",
        chrom="{chrom}",
        labels=" ".join(f"'{tool}'" for tool in TOOLS.values()),
        output_prefix=lambda wc, output: join(dirname(output[0]), wc.chrom),
    shell:
        """
        python {params.script} \
        -c {params.chrom} \
        -i {input.calls} \
        -l {params.labels} \
        -s {input.censat} \
        -d {input.segdup} \
        -o {params.output_prefix}
        """


rule sequence_context_all:
    input:
        expand(
            rules.plot_seq_context.output,
            chrom=CHROMS,
        ),
