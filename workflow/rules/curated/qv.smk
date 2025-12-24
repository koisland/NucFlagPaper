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


def cmd_contig_qv(wc, input):
    cmds = []
    # We don't know per contig QV filename since based on assembly filename.
    for file in input.merqury_qv:
        dname = dirname(file)
        fname, _ = splitext(basename(file))
        _, vnum = fname.rsplit("_", 1)
        qv_file = join(dname, f"{fname}.{vnum}.qv")
        cmds.append(f"<(cut -f 1,4 {qv_file})")

    return " ".join(cmds)


# QV correlation with NucFlag QV.
# * Should also show individual examples where QV underestimates (Show example in HSAT-1A in chr4_MATERNAL)
rule plot_nucflag_merqury_qv_plot:
    input:
        script="workflow/scripts/curated/qv_comparison.py",
        merqury_qv=expand(
            rules.qv_estimate_qv.output.qv,
            sm=[f"HG002_{v}" for v in reversed(ASSEMBLIES.keys())],
        ),
        nucflag_qv=expand(
            rules.calculate_nucflag_qv.output.qv, version=reversed(ASSEMBLIES.keys())
        ),
    output:
        plot=join(OUTPUT_DIR, "qv", "HG002_qv_{cond}.png"),
    conda:
        "../../envs/curated.yaml"
    params:
        # Merqury qv uses assembly name
        merqury_infiles=lambda wc, input: cmd_contig_qv(wc, input),
        nucflag_infiles=lambda wc, input: " ".join(
            f"<(cut -f 1,4 {file})" for file in input.nucflag_qv
        ),
        labels=" ".join(reversed(ASSEMBLIES.keys())),
        colors=" ".join(["blue", "green", "red"]),
        omit_acrocentrics=lambda wc: (
            "--omit-acrocentrics" if wc.cond == "omit_acro" else ""
        ),
    log:
        os.path.join(LOGS_DIR, "plot_nucflag_merqury_qv_plot_{cond}.log"),
    shell:
        """
        python {input.script} \
        -x {params.merqury_infiles} \
        -y {params.nucflag_infiles} \
        -l {params.labels} \
        -c {params.colors} \
        -o {output.plot} {params.omit_acrocentrics} 2> {log}
        """


rule qv_run_all:
    input:
        rules.qv_all.input,
        expand(rules.plot_nucflag_merqury_qv_plot.output, cond=["omit_acro", "default"]),
    default_target: True
