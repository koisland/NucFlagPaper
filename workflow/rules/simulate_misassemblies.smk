class DataUrlInfo(NamedTuple):
    url: str
    path: str


@dataclass
class DataSourceInfo:
    assembly: dict[str, DataUrlInfo] = field(default_factory=dict)
    hifi: dict[str, DataUrlInfo] = field(default_factory=dict)
    ont: dict[str, DataUrlInfo] = field(default_factory=dict)
    rm: dict[str, DataUrlInfo] = field(default_factory=dict)


DataInfo = tuple[list[str], list[str], list[str], defaultdict[str, DataSourceInfo]]
ALLOWED_DTYPES = {"assembly", "hifi", "ont", "rm"}
DATA_MANIFEST = config["simulate"]["manifest"]


def get_data_manifest(data_manifest: str, output_dir: str) -> DataInfo:
    samples = []
    dtypes = []
    url_hashes = []
    data = defaultdict(DataSourceInfo)
    with open(data_manifest) as fh:
        for line in fh:
            line = line.strip()
            sm, dtype, url = line.split()
            url_hash = hashlib.sha256(url.encode()).hexdigest()
            if dtype in ALLOWED_DTYPES:
                data_url_info = getattr(data[sm], dtype)
                base_fname = basename(url)
                data_url_info[url_hash] = DataUrlInfo(
                    url=url, path=join(output_dir, sm, dtype, base_fname)
                )
            else:
                raise ValueError(f"Invalid dtype ({dtype})")

            samples.append(sm)
            dtypes.append(dtype)
            url_hashes.append(url_hash)

    return samples, dtypes, url_hashes, data


def get_sample_misassembly_samples(
    data: defaultdict[str, DataSourceInfo]
) -> dict[str, dict[str, str | list[str]]]:
    samples = {}
    for sample, data_sources in data.items():
        for dtype, read_sources in [
            ("hifi", data_sources.hifi),
            ("ont_r10", data_sources.ont),
        ]:
            if not read_sources:
                continue

            data_sample_info = {}

            asm_info = next(iter(data_sources.assembly.values()))
            data_sample_info["asm_fa"] = asm_info.path
            data_sample_info["reads"] = [read.path for read in read_sources.values()]
            data_sample_info["config"] = config["config"][dtype]
            data_sample_info["flagger_config"] = config["config"][f"flagger_{dtype}"]
            data_sample_info["rm"] = [rm.path for rm in data_sources.rm.values()]
            data_sample_info["group_by"] = config["simulate"]["group_by"]
            samples[f"{sample}_{dtype}"] = data_sample_info

    return samples


def get_mtypes_w_seeds(
    number: int,
    lengths: list[int],
    mtypes: list[str],
) -> dict[int, dict[str, Any]]:
    return {
        # https://stackoverflow.com/a/67219726
        int.from_bytes(
            hashlib.sha256(f"{mtype}_{number}_{length}".encode()).digest()[:4], "little"
        ): {"mtype": mtype, "number": number, "length": length}
        for mtype in mtypes
        for length in lengths
    }


SAMPLES, DTYPES, URL_HASHES, DATA = get_data_manifest(DATA_MANIFEST, DATA_DIR)
NUM = 100
LEN = config.get("lengths", [1, 2, 10, 100, 1_000, 10_000, 50_000])
MTYPES = config.get("mtype", ["misjoin", "false_duplication", "inversion"])
DOWNSAMPLE_PERC = [0.50, 0.33]
ALL_MTYPES_SEEDED = get_mtypes_w_seeds(NUM, LEN, MTYPES)


wildcard_constraints:
    sm="|".join(SAMPLES),
    dtype="|".join(DTYPES),
    url_hash="|".join(URL_HASHES),


rule download_files:
    output:
        touch(join(DATA_DIR, "{sm}", "{dtype}", "download_{url_hash}.done")),
    log:
        join(LOGS_DIR, "{sm}", "{dtype}", "download_{url_hash}.log"),
    params:
        url=lambda wc: getattr(DATA[wc.sm], wc.dtype)[wc.url_hash].url,
        output_dir=lambda wc, output: dirname(output[0]),
    threads: 1
    conda:
        "../envs/download.yaml"
    shell:
        """
        if [[ "{params.url}" == s3* ]]; then
            aws s3 --no-sign-request cp --page-size 500 "{params.url}" {params.output_dir} 2> {log}
        else
            wget --no-verbose --no-check-certificate {params.url} -P {params.output_dir} 2> {log}
        fi
        """


rule download_all:
    input:
        expand(
            rules.download_files.output,
            zip,
            sm=SAMPLES,
            dtype=DTYPES,
            url_hash=URL_HASHES,
        ),


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
            "seeded_mtypes": ALL_MTYPES_SEEDED,
        }


use rule * from Misassemblies as sim_*
