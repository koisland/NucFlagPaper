

rule calculate_precision_recall:
    input:
        chkpt=[rules.ms_nucflag.input, rules.ms_all.input],
    output:
        summary=join(SUMMARY_OUTPUT_DIR, "nucflag_{sm}.tsv"),
        missed_calls_dir=directory(join(SUMMARY_OUTPUT_DIR, "nucflag_{sm}_missed")),
    params:
        script="workflow/scripts/metrics/calculate_precision_recall.py",
        fglob_test=lambda wc: f"{wc.sm}*_misassemblies.bed",
        fglob_truth=lambda wc: f"{wc.sm}/*.bed",
        dtype=get_dtype,
        output_dir=NUCFLAG_OUTPUT_DIR,
        sim_output_dir=MISASIM_OUTPUT_DIR,
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        python {params.script} \
        -a {params.output_dir} \
        --glob_test "{params.fglob_test}" \
        -b {params.sim_output_dir} \
        --glob_truth "{params.fglob_truth}" \
        -d {params.dtype} \
        --output_dir_missed_calls {output.missed_calls_dir} > {output.summary}
        """


use rule calculate_precision_recall as calculate_precision_recall_inspector with:
    input:
        chkpt=[rules.inspector.input, rules.ms_all.input],
    output:
        summary=join(SUMMARY_OUTPUT_DIR, "inspector_{sm}.tsv"),
        missed_calls_dir=directory(join(SUMMARY_OUTPUT_DIR, "inspector_{sm}_missed")),
    params:
        script="workflow/scripts/metrics/calculate_precision_recall.py",
        fglob_test=lambda wc: f"{wc.sm}*.bed",
        fglob_truth=lambda wc: f"{wc.sm}/*.bed",
        dtype=get_dtype,
        output_dir=join("results", "misasim", "inspector"),
        sim_output_dir=MISASIM_OUTPUT_DIR,


use rule calculate_precision_recall as calculate_precision_recall_flagger with:
    input:
        chkpt=rules.flagger.input,
    output:
        summary=join(SUMMARY_OUTPUT_DIR, "flagger_{sm}.tsv"),
        missed_calls_dir=directory(join(SUMMARY_OUTPUT_DIR, "flagger_{sm}_missed")),
    params:
        script="workflow/scripts/metrics/calculate_precision_recall.py",
        fglob_test=lambda wc: f"{wc.sm}*/*.bed",
        fglob_truth=lambda wc: f"{wc.sm}/*.bed",
        dtype=get_dtype,
        output_dir=join("results", "misasim", "flagger"),
        sim_output_dir=MISASIM_OUTPUT_DIR,


rule plot_precision_recall:
    input:
        summaries=[
            rules.calculate_precision_recall.output.summary,
            rules.calculate_precision_recall_inspector.output.summary,
            rules.calculate_precision_recall_flagger.output.summary,
        ],
    output:
        directory(join(SUMMARY_OUTPUT_DIR, "{sm}")),
    params:
        script="workflow/scripts/metrics/plot_precision_recall.py",
        output_prefix=lambda wc, output: f"{output[0]}",
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        python {params.script} -i {input.summaries} -o {params.output_prefix}/
        """


rule plot_precision_recall_xy:
    input:
        summaries=[
            rules.calculate_precision_recall.output.summary,
            rules.calculate_precision_recall_inspector.output.summary,
            rules.calculate_precision_recall_flagger.output.summary,
        ],
    output:
        directory(join(SUMMARY_OUTPUT_DIR, "{sm}_xy")),
    params:
        script="workflow/scripts/metrics/plot_precision_recall_xy.py",
        output_prefix=lambda wc, output: f"{output[0]}",
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        python {params.script} -i {input.summaries} -o {params.output_prefix}/
        """


rule plot_stats_f1_all:
    input:
        summaries=[
            rules.calculate_precision_recall.output.summary,
            rules.calculate_precision_recall_inspector.output.summary,
            rules.calculate_precision_recall_flagger.output.summary,
        ],
    output:
        tsv=join(SUMMARY_OUTPUT_DIR, "{sm}_summary_stats.tsv"),
        plot=join(SUMMARY_OUTPUT_DIR, "{sm}_summary_stats.png"),
    conda:
        "../../envs/tools.yaml"
    params:
        script="workflow/scripts/metrics/plot_summary_stats.py",
        colors=" ".join(TOOL_COLORS.values()),
        labels=" ".join(f"'{lbl}'" for lbl in TOOLS.values()),
    shell:
        """
        python {params.script} \
        -i {input.summaries} \
        -l {params.labels} \
        -c {params.colors} \
        -o {output.plot} > {output.tsv}
        """


rule generate_benchmarks:
    input:
        [
            rules.ms_nucflag.input,
            rules.flagger.input,
            rules.inspector.input,
        ],
    output:
        join(OUTPUT_DIR, "benchmark_plots", "rss.png"),
        join(OUTPUT_DIR, "benchmark_plots", "minutes.png"),
    params:
        benchmark_dir_args=" ".join(
            [
                f"{f} {join(BENCHMARK_DIR, tool)}"
                for (f, tool) in [
                    ("-i", "inspector"),
                    ("-f", "flagger"),
                    ("-i", "inspector"),
                ]
            ]
        ),
        script="workflow/scripts/metrics/aggregate_benchmarks.py",
        output_dir=lambda wc, output: dirname(output[0]),
    shell:
        """
        python {params.script} {params.benchmark_dir_args} -o {params.output_dir}
        """


rule benchmarks_all:
    input:
        expand(rules.calculate_precision_recall.output, sm=SAMPLES),
        expand(rules.calculate_precision_recall_inspector.output, sm=SAMPLES),
        expand(rules.calculate_precision_recall_flagger.output, sm=SAMPLES),
        expand(rules.plot_precision_recall.output, sm=SAMPLES),
        expand(rules.plot_precision_recall_xy.output, sm=SAMPLES),
        expand(rules.plot_stats_f1_all.output, sm=SAMPLES),
        rules.generate_benchmarks.output,
