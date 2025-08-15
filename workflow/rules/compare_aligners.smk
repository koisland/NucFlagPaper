

module AlignerComparison:
    snakefile:
        "aligners/Snakefile"
    config:
        {
            "samples": get_sample_misassembly_samples({"HG002": DATA["HG002"]}),
            "output_dir": join(config["output_dir"], "aligners"),
            "logs_dir": join(config["logs_dir"], "aligners"),
            "benchmarks_dir": join(config["benchmarks_dir"], "aligners"),
            "data": [v.path for v in DATA["HG002"].hifi.values()],
            "misasim_output_dir": join(config["output_dir"], "misasim"),
            "seeded_mtypes": ALL_MTYPES_SEEDED,
        }


use rule * from AlignerComparison as aln_cmp_*
