"""
{
    Literal["output_dir"]: str,
    Literal["logs_dir"]: str,
    Literal["benchmarks_dir"]: str,
    Literal["samples"]: {
        str: {
            Literal["asm_fa"]: str,
            Literal["bed"]: str,
            Literal["reads"]: list[str],
            Literal["preset"]: str
        }
    }
}
"""

module Misassemblies:
    snakefile:
        "misasim/Snakefile"
    config:
        {
            "output_dir": join(config["output_dir"], "misasim"),
            "logs_dir": join(config["logs_dir"], "misasim"),
            "benchmarks_dir": join(config["benchmarks_dir"], "misasim"),
            "samples": get_sample_misassembly_samples(DATA),
            "downsample_perc": DOWNSAMPLE_PERC,
            "seeded_mtypes": ALL_MTYPES_SEEDED
        }


use rule * from Misassemblies as sim_*
