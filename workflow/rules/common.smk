import re
import os
import hashlib
import json
from os.path import join, dirname, basename, splitext
from typing import NamedTuple, Any
from collections import defaultdict
from dataclasses import dataclass, field


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


LOGS_DIR = config.get("logs_dir", "logs")
BMKS_DIR = config.get("benchmarks_dir", "benchmarks")
DATA_DIR = config.get("data_dir", "data")
DATA_MANIFEST = config.get("data_manifest")
RGX_DTYPE = re.compile(r"^(.*?)_(hifi|ont_r9|ont_r10)")


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
            data_sample_info["group_by"] = config["group_by"]
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


def get_dtype(wc) -> str:
    mtch = RGX_DTYPE.search(wc.sm)
    if not mtch:
        raise ValueError(f"{wc.sm} doesn't match regex ({RGX_DTYPE})")
    sm, dtype = mtch.groups()
    return dtype
