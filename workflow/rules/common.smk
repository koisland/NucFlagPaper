import re
import os
import hashlib
import json
from os.path import join, dirname, basename, splitext
from typing import NamedTuple, Any
from collections import defaultdict
from dataclasses import dataclass, field


LOGS_DIR = config.get("logs_dir", "logs")
BMKS_DIR = config.get("benchmarks_dir", "benchmarks")
DATA_DIR = config.get("data_dir", "data")
RGX_DTYPE = re.compile(r"^(.*?)_(hifi|ont_r9|ont_r10)")


def get_dtype(wc) -> str:
    mtch = RGX_DTYPE.search(wc.sm)
    if not mtch:
        raise ValueError(f"{wc.sm} doesn't match regex ({RGX_DTYPE})")
    sm, dtype = mtch.groups()
    return dtype
