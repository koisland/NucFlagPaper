import yaml


with open(config["assembly_config"], "rb") as fh:
    assembly_config = yaml.safe_load(fh)


module CompareAssembly:
    snakefile:
        "compare_assembly/Snakefile"
    config:
        {
            "samples": assembly_config["samples"],
            "output_dir": join(config["output_dir"], "assembly"),
            "logs_dir": join(config["logs_dir"], "assembly"),
            "benchmarks_dir": join(config["benchmarks_dir"], "assembly"),
        }


use rule * from CompareAssembly as cmp_asm_*
