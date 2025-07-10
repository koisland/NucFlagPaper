
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
    default_target: True
