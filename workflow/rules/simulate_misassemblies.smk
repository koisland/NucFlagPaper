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


def get_sample_misassembly_samples(
    data: defaultdict[str, DataSourceInfo]
) -> dict[str, dict[str, str | list[str]]]:
    samples = {}
    for sample, data_sources in data.items():
        for dtype, read_sources in [
            ("hifi", data_sources.hifi),
            ("ont", data_sources.ont),
        ]:
            data_sample_info = {}

            asm_info = next(iter(data_sources.assembly.values()))
            data_sample_info["asm_fa"] = asm_info.path
            data_sample_info["reads"] = [read.path for read in read_sources.values()]
            data_sample_info["bed"] = config["regions"][sample]
            data_sample_info["preset"] = dtype

            samples[sample] = data_sample_info

    return samples


module Misassemblies:
    snakefile:
        "misasim/Snakefile"
    config:
        {
            "output_dir": config["output_dir"],
            "logs_dir": config["logs_dir"],
            "benchmarks_dir": config["benchmarks_dir"],
            "samples": get_sample_misassembly_samples(DATA),
        }


use rule * from Misassemblies as sim_*
