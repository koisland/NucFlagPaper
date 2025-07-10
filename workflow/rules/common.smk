import hashlib
from os.path import join, dirname, basename
from typing import NamedTuple
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


DataInfo = tuple[list[str], list[str], list[str], defaultdict[str, DataSourceInfo]]
ALLOWED_DTYPES = {"assembly", "hifi", "ont"}


LOGS_DIR = config.get("logs_dir", "logs")
BMKS_DIR = config.get("benchmarks_dir", "benchmarks")
DATA_DIR = config.get("data_dir", "data")
DATA_MANIFEST = config["data_manifest"]


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


SAMPLES, DTYPES, URL_HASHES, DATA = get_data_manifest(DATA_MANIFEST, DATA_DIR)


wildcard_constraints:
    sm="|".join(SAMPLES),
    dtype="|".join(DTYPES),
    url_hash="|".join(URL_HASHES),
