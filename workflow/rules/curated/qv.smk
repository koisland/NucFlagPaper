QV_CONFIG = {
    "output_dir": join(OUTPUT_DIR, "qv"),
    "log_dir": join(LOG_DIR, "qv"),
    "benchmark_dir": join(BENCHMARK_DIR, "qv"),
    "samples": {
        f"HG002_{version}": {
            "asm": expand(rules.download_curated_asm.output.fa, version=version),
            "reads": rules.download_element.output,
        }
        for version in ASSEMBLIES.keys()
    },
    "threads": 12,
    "mem": "20GB",
    "kmer_size": 31,
}


module QV:
    snakefile:
        "Snakemake-merqury/workflow/Snakefile"
    config:
        QV_CONFIG


use rule * from QV as qv_*


rule calculate_nucflag_qv:
    input:
        calls=expand(
            rules.versioned_check_asm_nucflag.output.misassemblies,
            sm="HG002_{version}",
        ),
    output:
        qv=join(OUTPUT_DIR, "qv", "HG002_{version}_nucflag_qv.bed"),
    conda:
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    log:
        os.path.join(LOGS_DIR, "run_nucflag_element_HG002_{version}.log"),
    shell:
        """
        nucflag qv -i {input.calls} > {output.qv} 2> {log}
        """


# TODO: QV correlation with NucFlag QV.
# * Should also show that QV underestimates (Show example in HSAT-1A in chr4_MATERNAL)


rule plot_nucflag_merqury_qv_plot:
    input:
        script="workflow/scripts/curated/qv_comparison.py",
        merqury_qv=expand(
            rules.qv_estimate_qv.output.qv,
            sm=[f"HG002_{v}" for v in ASSEMBLIES.keys()],
        ),
        nucflag_qv=expand(
            rules.calculate_nucflag_qv.output.qv, version=ASSEMBLIES.keys()
        ),
    output:
        plot=join(OUTPUT_DIR, "qv", "HG002_qv.png"),
    conda:
        "../../envs/curated.yaml"
    params:
        merqury_infiles=lambda wc, input: " ".join(
            f"<(cut -f 1,4 {file})" for file in input.merqury_qv
        ),
        nucflag_infiles=lambda wc, input: " ".join(
            f"<(cut -f 1,4 {file})" for file in input.nucflag_qv
        ),
        labels=" ".join(ASSEMBLIES.keys()),
        colors=" ".join(["red", "green", "blue"]),
    log:
        os.path.join(LOGS_DIR, "plot_nucflag_merqury_qv_plot.log"),
    shell:
        """
        python {input.script} \
        -x {params.merqury_infiles} \
        -y {params.nucflag_infiles} \
        -l {params.labels} \
        -c {params.colors} \
        -o {output.plot} 2> {log}
        """


use rule qv_all from QV as qv_all with:
    input:
        rules.qv_all.input,
        rules.plot_nucflag_merqury_qv_plot.output,
    default_target: True
