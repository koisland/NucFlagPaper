
rule convert_vcf_to_bed:
    input:
        rules.download_vcf.output,
    output:
        os.path.join(OUTPUT_DIR, "data", "{version}_truth.bed"),
    shell:
        """
        zcat -f {input} | grep -v "#" | awk -v OFS="\\t" '{{
            print $1, $2, $2 + length($4), $4"-"$5, 0, $3, $2, $2 + length($4), "0,0,0"
        }}' > {output}
        """


rule calculate_precision_recall:
    input:
        test_bed=get_bed_files,
        truth_bed=rules.convert_vcf_to_bed.output,
    output:
        summary=join(OUTPUT_DIR, "summary", "{tool}_{version}.tsv"),
        missed_calls_dir=directory(
            join(OUTPUT_DIR, "summary", "{tool}_{version}_missed")
        ),
    params:
        script="workflow/scripts/metrics/calculate_precision_recall_curated.py",
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        python {params.script} \
        -a {input.test_bed} \
        -b {input.truth_bed} \
        --output_dir_missed_calls {output.missed_calls_dir} > {output.summary}
        """


# Create venn-diagram shared between callers
# Shows that NucFlag picks up a lot more than other callers.
rule plot_call_venn:
    input:
        nucflag=expand(
            rules.versioned_check_asm_nucflag.output.misassemblies,
            sm="HG002_{version}",
        )[0],
        flagger=expand(rules.versioned_run_flagger.output, sm="HG002_{version}")[0],
        inspector=expand(rules.versioned_merge_calls.output, sm="HG002_{version}")[0],
        deepvariant=expand(
            rules.versioned_filter_vcf_calls.output, sm="HG002_{version}"
        )[0],
    output:
        plot=join(OUTPUT_DIR, "call_ovl", "{version}_venn.png"),
    params:
        script="workflow/scripts/curated/call_venn_diagram.py",
        output_prefix=lambda wc, output: join(dirname(output[0]), wc.version),
        json_inputs=lambda wc, input: json.dumps(dict(input)),
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        python {params.script} -i '{params.json_inputs}' -o {params.output_prefix}
        """


rule plot_stats_f1_all:
    input:
        summaries=[
            expand(
                rules.calculate_precision_recall.output.summary,
                tool=tool,
                version=version,
            )
            for tool in TOOLS.keys()
            for version in ASSEMBLIES.keys()
        ],
    output:
        tsv=join(OUTPUT_DIR, "summary", "all_summary_stats.tsv"),
        plot=join(OUTPUT_DIR, "summary", "all_summary_stats.png"),
    conda:
        "../../envs/tools.yaml"
    params:
        script="workflow/scripts/metrics/plot_summary_stats.py",
        colors=" ".join(
            TOOL_COLORS[tool] for tool in TOOLS.keys() for _ in ASSEMBLIES.keys()
        ),
        labels=" ".join(
            f"'{TOOLS[tool]} ({version})'"
            for tool in TOOLS.keys()
            for version in ASSEMBLIES.keys()
        ),
    shell:
        """
        python {params.script} \
        -i {input.summaries} \
        -l {params.labels} \
        -c {params.colors} \
        -o {output.plot} > {output.tsv}
        """


rule benchmarks_all:
    input:
        expand(rules.convert_vcf_to_bed.output, version=ASSEMBLIES.keys()),
        expand(rules.plot_call_venn.output, version=ASSEMBLIES.keys()),
        expand(
            rules.calculate_precision_recall.output,
            version=ASSEMBLIES.keys(),
            tool=TOOLS.keys(),
        ),
        rules.plot_stats_f1_all.output,
