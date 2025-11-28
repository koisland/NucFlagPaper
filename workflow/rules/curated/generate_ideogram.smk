# Download data source
rule download_annotation_data:
    output:
        segdups=join(IDEOGRAM_OUTPUT_DIR, "HG002.SDs.010624.45col.bb"),
        censat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.0.1.cenSatv2.0.noheader.bb"),
        cytoband=join(IDEOGRAM_OUTPUT_DIR, "cytoBand.hg002v1.0.bb"),
    params:
        output_dir=IDEOGRAM_OUTPUT_DIR,
        urls=" ".join(ANNOTATIONS_V101),
    shell:
        """
        for url in {params.urls}; do
            fname=$(basename "${{url}}")
            new_path="{params.output_dir}/${{fname}}"
            wget "${{url}}" -O "${{new_path}}"
        done
        """


rule convert_to_bed:
    input:
        segdups=rules.download_annotation_data.output.segdups,
        censat=rules.download_annotation_data.output.censat,
        cytoband=rules.download_annotation_data.output.cytoband,
    output:
        segdups=join(IDEOGRAM_OUTPUT_DIR, "HG002.SDs.010624.45col.bed"),
        censat=join(IDEOGRAM_OUTPUT_DIR, "hg002v1.0.1.cenSatv2.0.noheader.bed"),
        cytoband=join(IDEOGRAM_OUTPUT_DIR, "cytoBand.hg002v1.0.bed"),
    shell:
        """
        bigbedtobed {input.segdups} {output.segdups}
        bigbedtobed {input.censat} {output.censat}
        bigbedtobed {input.cytoband} {output.cytoband}
        """


rule format_cytoband:
    input:
        censat=rules.convert_to_bed.output.censat,
        cytoband=rules.convert_to_bed.output.cytoband,
        fai=expand(rules.download_curated_asm.output.fai, version=f"v1.0.1")[0],
    output:
        cytoband=join(IDEOGRAM_OUTPUT_DIR, "cytoBand.hg002v1.0.formatted.bed"),
    params:
        bp_merge=5_000_000,
    shell:
        """
        cat \
            <(grep active_hor {input.censat} | \
            bedtools merge -i - -d {params.bp_merge} | \
            awk -v OFS="\\t" '{{
                midpt=int((($3-$2)/2) + $2);
                print $1, $2, midpt, ".", "acen";
                print $1, midpt, $3, ".", "acen"
            }}') \
            {input.cytoband} | \
        sort -k1,1 -k2,2n > {output.cytoband}.tmp

        bedtools subtract \
            -a <(grep -v <(printf "chrEBV\\nchrM\\n") {input.fai} | awk -v OFS="\\t" '{{ print $1, 0, $2 }}') \
            -b {output.cytoband}.tmp | \
        awk -v OFS="\\t" '{{ print $0, ".", "gneg" }}' >> {output.cytoband}.tmp
        sort -k1,1 -k2,2n {output.cytoband}.tmp > {output.cytoband}
        rm -f {output.cytoband}.tmp
        """


rule plot_ideogram:
    input:
        script="workflow/scripts/metrics/plot_curated_ideogram.py",
        cytobands=rules.format_cytoband.output,
        truth=expand(rules.convert_vcf_to_bed.output, version=f"v1.0.1"),
        nucflag_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="nucflag",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        inspector_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="inspector",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        flagger_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="flagger",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        fai=expand(rules.download_curated_asm.output.fai, version=f"v1.0.1")[0],
        segdups=rules.convert_to_bed.output.segdups,
        censat=rules.convert_to_bed.output.censat,
    output:
        plots=[
            join(IDEOGRAM_OUTPUT_DIR, "ideogram.png"),
            join(IDEOGRAM_OUTPUT_DIR, "ideogram_overview.png"),
            join(IDEOGRAM_OUTPUT_DIR, "ideogram.pdf"),
        ],
    params:
        output_prefix=lambda wc, output: os.path.splitext(output.plots[0])[0],
    shell:
        """
        python {input.script} \
        --fai {input.fai} \
        --cytobands {input.cytobands} \
        --truth {input.truth} \
        --nucflag {input.nucflag_calls} \
        --inspector {input.inspector_calls} \
        --flagger {input.flagger_calls} \
        --segdups {input.segdups} \
        --censat {input.censat} \
        -o {params.output_prefix}
        """


rule plot_ideogram_chrom:
    input:
        script="workflow/scripts/metrics/plot_curated_ideogram_chrom.py",
        cytobands=rules.format_cytoband.output,
        truth=expand(rules.convert_vcf_to_bed.output, version=f"v1.0.1"),
        nucflag_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="nucflag",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        inspector_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="inspector",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        flagger_calls=join(
            expand(
                rules.calculate_precision_recall.output.missed_calls_dir,
                tool="flagger",
                version=f"v1.0.1",
            )[0],
            "missed_calls.tsv",
        ),
        segdups=rules.convert_to_bed.output.segdups,
        censat=rules.convert_to_bed.output.censat,
    output:
        plot=join(IDEOGRAM_OUTPUT_DIR, "{chrom}.png"),
    shell:
        """
        python {input.script} \
        --cytobands {input.cytobands} \
        --truth {input.truth} \
        --nucflag {input.nucflag_calls} \
        --inspector {input.inspector_calls} \
        --flagger {input.flagger_calls} \
        --segdups {input.segdups} \
        --censat {input.censat} \
        -o {output.plot} \
        -c {wildcards.chrom}
        """


rule generate_ideogram:
    input:
        rules.plot_ideogram.output,
        expand(
            rules.plot_ideogram_chrom.output,
            chrom=["chrY_PATERNAL", "chr21_PATERNAL", "chr4_PATERNAL"],
        ),
    default_target: True
