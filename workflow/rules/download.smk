class DataUrlInfo(NamedTuple):
    url: str
    path: str
    # Command to run on data.
    cmd: str | None


@dataclass
class DataSourceInfo:
    assembly: dict[str, DataUrlInfo] = field(default_factory=dict)
    hifi: dict[str, DataUrlInfo] = field(default_factory=dict)
    ont: dict[str, DataUrlInfo] = field(default_factory=dict)
    rm: dict[str, DataUrlInfo] = field(default_factory=dict)
    vcf: dict[str, DataUrlInfo] = field(default_factory=dict)
    cytobands: dict[str, DataUrlInfo] = field(default_factory=dict)

    def get_first(self, dtype: str) -> DataUrlInfo:
        return next(iter(getattr(self, dtype).values()))

class DataInfo(NamedTuple):
    samples: list[str]
    dtypes: list[str]
    hashes: list[str]
    data: DataSourceInfo

ALLOWED_DTYPES = {"assembly", "hifi", "ont", "rm", "vcf", "cytobands"}
DATA_MANIFEST = config["simulate"]["manifest"]


def get_data_manifest(data_manifest: str, output_dir: str) -> DataInfo:
    samples = []
    dtypes = []
    url_hashes = []
    data = defaultdict(DataSourceInfo)
    with open(data_manifest) as fh:
        for line in fh:
            line = line.strip()
            try:
                sm, dtype, url, cmd = line.split("\t")
            except ValueError:
                sm, dtype, url = line.split("\t")
                cmd = None

            url_hash = hashlib.sha256(url.encode()).hexdigest()

            if dtype in ALLOWED_DTYPES:
                data_url_info = getattr(data[sm], dtype)
                base_fname = basename(url)
                data_url_info[url_hash] = DataUrlInfo(
                    url=url,
                    path=join(output_dir, sm, dtype, base_fname),
                    cmd=cmd if cmd else None
                )
            else:
                raise ValueError(f"Invalid dtype ({dtype})")

            samples.append(sm)
            dtypes.append(dtype)
            url_hashes.append(url_hash)

    return DataInfo(samples, dtypes, url_hashes, data)

# Get data
DATA_INFO = get_data_manifest(DATA_MANIFEST, DATA_DIR)
SAMPLES, DTYPES, URL_HASHES, DATA = DATA_INFO

# Generate list of expected files after download.
EXPECTED_FILES = []
for sm, dtype, hsh in zip(SAMPLES, DTYPES, URL_HASHES, strict=True):
    path = getattr(DATA[sm], dtype)[hsh].path
    EXPECTED_FILES.append(path)        
    if dtype == "assembly":
        EXPECTED_FILES.append(f"{path}.fai")   


rule download_files:
    output:
        touch(join(DATA_DIR, "{sm}", "{dtype}", "download_{url_hash}.done")),
    log:
        join(LOGS_DIR, "{sm}", "{dtype}", "download_{url_hash}.log"),
    params:
        url=lambda wc: getattr(DATA[wc.sm], wc.dtype)[wc.url_hash].url,
        path=lambda wc: getattr(DATA[wc.sm], wc.dtype)[wc.url_hash].path,
        cmd=lambda wc: getattr(DATA[wc.sm], wc.dtype)[wc.url_hash].cmd,
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

        # If SRA, convert to fastx
        if [[ "{params.path}" == "*sra*" ]]; then
            fastq-dump {params.path} --stdout > {params.path}.tmp
            rm -f {params.path} && mv {params.path}.tmp {params.path}
        fi

        # Apply command
        if [[ "{params.cmd}" != "" ]]; then
            {params.cmd} {params.path} > {params.path}.tmp
            rm -f {params.path} && mv {params.path}.tmp {params.path}
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
        )
    output:
        # Files that should exist at this point.
        EXPECTED_FILES
