

rule aggregate_benchmarks:
    input:
        chkpt=rules.plot_err_comparison.output,
        bmk_dir="benchmarks/hprc",
    output:
        plot=join(OUTPUT_DIR, "benchmarks", "hprc_runtime_metrics.png"),
    params:
        script="workflow/scripts/compare_hprc/aggregate_benchmarks.py",
        output_prefix=lambda wc, output: splitext(output.plot)[0],
    shell:
        """
        python {params.script} -i {input.bmk_dir} -o {params.output_prefix}
        """
