import yaml


with open(config["assembly_config"]["manifest"], "rb") as fh:
    assembly_config = yaml.safe_load(fh)
    sample_config = assembly_config["samples"]
    for sm, cfg in sample_config.items():
        data = cfg["data"]
        data["ont"]["path"] = [data.path for data in DATA["CHM13"].ont.values()]
        data["hifi"]["path"] = [data.path for data in DATA["CHM13"].hifi.values()]


module CompareAssembly:
    snakefile:
        "compare_assembly/Snakefile"
    config:
        {
            "samples": assembly_config["samples"],
            "output_dir": join(config["output_dir"], "assembly"),
            "log_dir": join(config["logs_dir"], "assembly"),
            "benchmark_dir": join(config["benchmarks_dir"], "assembly"),
            "nucflag_config": {
                "hifi": config["assembly_config"]["nucflag_hifi"],
                "ont": config["assembly_config"]["nucflag_ont"],
            },
            "asm_chm13": DATA["CHM13"].get_first("assembly").path,
        }


use rule * from CompareAssembly as cmp_asm_*
