"""
Rename the assemblies. We have two verkko versions so need a prefix.
"""


rule rename_assemblies:
    input:
        asm=lambda wc: (
            expand(
                rules.denovo_verkko_output.output, asm=wc.sm.split("_")[1], sm=wc.sm
            )[0]
            if wc.sm.split("_")[1] == "verkko"
            else expand(
                rules.denovo_hifiasm_output.output, asm=wc.sm.split("_")[1], sm=wc.sm
            )[0]
        ),
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
        join - <(cut -f1,2 {input.fai}) | \
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
            sm=config["samples"].keys(),
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
        expand(rules.rename_assemblies.output.asm, sm=sm)
        for sm in config["samples"].keys()
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


use rule all from VariantGraph as vg_all with:
    input:
        rules.vg_all.input,
        expand(rules.intersect_bubbles_w_calls.output, sm=config["samples"].keys()),
    default_target: True
