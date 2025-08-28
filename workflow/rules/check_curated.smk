module HG002_curated:
    snakefile:
        "curated/Snakefile"
    config:
        {
            "output_dir": join(config["output_dir"], "curated"),
            "logs_dir": join(config["logs_dir"], "curated"),
            "benchmarks_dir": join(config["benchmarks_dir"], "curated"),
            "hifi": [v.path for v in DATA["HG002"].hifi.values()],
            "config": config["config"]["hifi_curated"],
            "alpha": config["config"]["flagger_hifi"],
        }


use rule * from HG002_curated as hg002_curated_*
