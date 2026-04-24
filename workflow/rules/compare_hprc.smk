module CompareHPRCReleases:
    snakefile:
        "compare_hprc/Snakefile"
    config:
        {
            "output_dir": join(config["output_dir"], "hprc"),
            "logs_dir": join(config["logs_dir"], "hprc"),
            "benchmarks_dir": join(config["benchmarks_dir"], "hprc"),
            "manifest": config["hprc"]["manifest"],
            "annotations_segdups": {
                "R1": config["hprc"]["annotations_segdups"]["R1"],
                "R2": config["hprc"]["annotations_segdups"]["R2"],
            },
        }


use rule * from CompareHPRCReleases as hprc_*
