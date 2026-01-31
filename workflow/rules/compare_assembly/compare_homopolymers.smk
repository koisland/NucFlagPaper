rule find_all_lcr_regions:
    input:
        ref=lambda wc: get_assembly(wc.sm),
    output:
        bed=join(OUTPUT_DIR, "homopolymers", "{sm}_lcr.bed"),
    conda:
        "../../envs/curated.yaml"
    log:
        join(LOG_DIR, "homopolymers", "find_all_lcr_regions_{sm}.log"),
    shell:
        """
        longdust {input.ref} > {output.bed} 2> {log}
        """


# What distribution of nucflag calls have abnormal signals
# Where do these errors fall in across all homopolymers?
rule label_lcr_regions_and_pos:
    input:
        script="workflow/scripts/curated/label_homopolymers.py",
        ref=lambda wc: get_assembly(wc.sm),
        bed=rules.find_all_lcr_regions.output,
        calls=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm="{sm}_hifi",
        ),
    output:
        bed=join(OUTPUT_DIR, "homopolymers", "{sm}_homopolymers.bed.gz"),
        bins=join(OUTPUT_DIR, "homopolymers", "{sm}_bin_cnts.tsv"),
        plots=[
            join(OUTPUT_DIR, "homopolymers", "{sm}_ovl_stacked.png"),
            join(OUTPUT_DIR, "homopolymers", "{sm}_ovl_split.png"),
            join(OUTPUT_DIR, "homopolymers", "{sm}_all_stacked.png"),
            join(OUTPUT_DIR, "homopolymers", "{sm}_all_split.png"),
        ],
    params:
        plot_prefix=lambda wc, output: join(dirname(output.plots[0]), wc.sm),
    conda:
        "../../envs/curated.yaml"
    log:
        join(LOG_DIR, "homopolymers", "label_lcr_regions_and_pos_{sm}.log"),
    shell:
        """
        python {input.script} \
        -i {input.bed} \
        -r {input.ref} \
        -c {input.calls} \
        -b {output.bins} \
        -p {params.plot_prefix} | bgzip > {output.bed}
        """


rule overlap_cmp_dist:
    input:
        script="workflow/scripts/compare_assembly/cmp_homopolymer_dist.py",
        bins=expand(
            rules.label_lcr_regions_and_pos.output.bins, sm=config["samples"].keys()
        ),
    output:
        plot=join(OUTPUT_DIR, "homopolymers", "dist_homopolymers_w_signif.png"),
    conda:
        "../../envs/curated.yaml"
    params:
        labels=" ".join(config["samples"].keys()),
    shell:
        """
        python {input.script} \
        --bins {input.bins} \
        --labels {params.labels} \
        -o {output.plot}
        """


rule cmp_homopolymers_all:
    input:
        expand(rules.label_lcr_regions_and_pos.output, sm=config["samples"].keys()),
        rules.overlap_cmp_dist.output,
