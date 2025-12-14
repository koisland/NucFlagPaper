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
}


module QV:
    snakefile:
        "Snakemake-merqury/workflow/Snakefile"
    config:
        QV_CONFIG


use rule * from QV as qv_*
