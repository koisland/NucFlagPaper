VG_SAMPLES = ["CHM13v2.0", *config["samples"].keys()]


"""
Rename the assemblies. We have two verkko versions so need a prefix.
"""


rule rename_assemblies:
    input:
        asm=lambda wc: get_assembly(wc.sm),
    output:
        asm=join(OUTPUT_DIR, "vg", "{sm}.fa.gz"),
        asm_fai=join(OUTPUT_DIR, "vg", "{sm}.fa.gz.fai"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        seqkit replace -p '^' -r '{wildcards.sm}_' {input.asm} | bgzip > {output.asm}
        samtools faidx {output.asm}
        """


"""
Generate a filtered callset from both hifi.
Rename bed to match assembly.
Filter small regions.
"""


rule generate_filtered_calls:
    input:
        beds=lambda wc: expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm=f"{wc.sm}_hifi",
        ),
        fai=rules.rename_assemblies.output.asm_fai,
    output:
        join(OUTPUT_DIR, "vg", "{sm}_filtered_calls.bed"),
    params:
        bp_slop=5_000,
        bp_ignore_ends=100_000,
        min_ctg_length=1_000_000,
        rgx_call_filter="insertion|deletion|false_dup|misjoin",
    # conda:
    #     "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    conda:
        "../../envs/compare_assembly.yaml"
    log:
        join(LOG_DIR, "vg", "{sm}_consensus.log"),
    shell:
        """
        {{ grep -hP "{params.rgx_call_filter}" {input.beds} | \
        awk -v OFS="\\t" '{{ $1="{wildcards.sm}_"$1; print }}' | \
        join - <(cut -f1,2 {input.fai} | sort -k1,1 -k2,2n) | \
        awk -v OFS="\\t" '{{
            ctg_len=$(NF);
            len=$3-$2;
            ovl_st=$3 < {params.bp_ignore_ends}
            ovl_end=$2 > (ctg_len - {params.bp_ignore_ends}) && $3 < ctg_len
            if (ctg_len < {params.min_ctg_length} || ovl_st || ovl_end) {{ next; }};
            print $1, $2, $3, $4
        }}' | \
        bedtools slop -i - -b {params.bp_slop} -g {input.fai} ;}} > {output} 2> {log}
        """


checkpoint merge_beds:
    input:
        beds=expand(
            rules.generate_filtered_calls.output,
            sm=VG_SAMPLES,
        ),
    output:
        calls_merged=join(OUTPUT_DIR, "vg", "filtered_calls_merged.bed"),
    conda:
        "../../envs/compare_assembly.yaml"
    params:
        merge_by=5_000,
        min_rgn_length=5_000,
    shell:
        """
        bedtools merge -i <(cat {input.beds}) -d {params.merge_by} -c 4 -o collapse | \
        awk '$3-$2 > {params.min_rgn_length}' > {output.calls_merged}
        """


ASM_VAR_GRAPH_CONFIG = {
    "assemblies": [
        expand(rules.rename_assemblies.output.asm, sm=sm) for sm in VG_SAMPLES
    ],
    "regions": rules.merge_beds.output.calls_merged,
    "output_dir": join(OUTPUT_DIR, "vg"),
    "log_dir": join(LOG_DIR, "vg"),
    "benchmark_dir": join(BENCHMARK_DIR, "vg"),
    "threads_mm2": 12,
    "mem_mm2": "100GB",
    "threads_mg": 12,
    "mem_mg": "30GB",
}


# Align hifiasm to verkko.
module VariantGraph:
    snakefile:
        "Snakemake-asm-evaluation-vg/workflow/Snakefile"
    config:
        ASM_VAR_GRAPH_CONFIG


use rule * from VariantGraph as vg_*


"""
Intersect bubbles in variation graph with calls to pinpoint regions to show.
Generally low consensus regions or errors unique to an assembly.
Visualize with Bandage and IGV
"""


# TODO: Ratio between small and large allele in bubbles to pick up larger errors. Also plot.
rule intersect_bubbles_w_calls:
    input:
        bed=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm="{sm}_hifi",
        ),
        bubbles=rules.vg_generate_variation_graph_output.output,
    output:
        bed_intersect=join(OUTPUT_DIR, "vg", "{sm}_bubble_call_intersection.bed"),
    conda:
        "../../envs/compare_assembly.yaml"
    shell:
        """
        # Add sample name to chrom col based on filename to match bubbles bed.
        bedtools intersect \
        -a <(awk -v OFS="\\t" 'FNR > 1 {{
            $1="{wildcards.sm}_"$1;
            print
        }}' {input.bed} | sort -k1,1 -k2,2n) \
        -b {input.bubbles} \
        -wa -wb | \
        awk -v OFS="\\t" '{{ sub("{wildcards.sm}_", "", $1); print }}' > {output.bed_intersect}
        """


# Liftover relative to CHM13v2.0 across all genomic intervals (excluding chrY).
rule query_chm13_impg:
    input:
        paf=rules.vg_ava_asm.output,
        fai=rules.download_chm13.output.fai,
    output:
        bed=join(OUTPUT_DIR, "vg", "ava_asm_liftover_chm13.bed.gz"),
    params:
        chm13_prefix=f"{REF_SM}_",
        window=100_000,
    conda:
        "Snakemake-asm-evaluation-vg/workflow/envs/tools.yaml"
    shell:
        """
        # Make windows for chm13v2.0
        # Query all regions from ava alignment.
        impg query \
        -b <(bedtools makewindows -b <(awk -v OFS="\\t" '{{ print "{params.chm13_prefix}"$1, 0, $2 }}' {input.fai} | grep -v chrY) -w {params.window}) \
        -p {input.paf} |
        gzip > {output.bed}
        """


rule plot_heatmap_chm13_impg:
    input:
        bed=rules.query_chm13_impg.output,
        calls=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm=[f"{sm}_hifi" for sm in VG_SAMPLES],
        ),
        annots=[
            rules.download_annotations.output.censat,
            rules.download_annotations.output.segdups,
        ],
    output:
        plot_dir=directory(join(OUTPUT_DIR, "vg", "plot_ava_asm_liftover_chm13")),
    params:
        script="workflow/scripts/compare_assembly/plot_impg_query_heatmap.py",
        annots=lambda wc, input: " ".join(
            [
                f"""<(awk -v OFS='\\t' '{{ $1="{REF_SM}_"$1; print }}' {annot})"""
                for annot in input.annots
            ]
        ),
        annot_labels=["'Satellite structure'", "'Segmental duplications'"],
        labels=" ".join([REF_SM, *config["samples"].keys()]),
        colors=" ".join(["red", *[ASM_COLORS[sm] for sm in config["samples"].keys()]]),
    shell:
        """
        python {params.script} \
        -i {input.bed} \
        -b {input.calls} \
        -l {params.labels} \
        -c {params.colors} \
        -a {params.annots} \
        -al {params.annot_labels} \
        -o {output.plot_dir}
        """


use rule all from VariantGraph as vg_all with:
    input:
        rules.vg_all.input,
        rules.plot_heatmap_chm13_impg.output,
        expand(rules.intersect_bubbles_w_calls.output, sm=VG_SAMPLES),
    default_target: True
