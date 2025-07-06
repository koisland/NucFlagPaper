from os.path import join


THREADS = 4
GLOB_DATA_MANIFEST = "{exp}/manifest.tsv"
EXPERIMENTS = glob_wildcards(GLOB_DATA_MANIFEST).exp


wildcard_constraints:
    exp="|".join(EXPERIMENTS)


rule download_manifest:
    input:
        script="download.sh",
        manifest=GLOB_DATA_MANIFEST,
    output:
        touch(join("{exp}", "download_data.done"))
    log:
        "logs/download_manifest_{exp}.log"
    threads: THREADS
    conda:
        "env.yaml"
    shell:
        """
        awk -v OFS='\\t' '{{ print $2, $3 }}' {input.manifest} | \
        parallel -j {threads} --colsep "\\t" 'bash {input.script} {{2}} {wildcards.exp}/{{1}}' 2> {log}
        """


rule download_all:
    input:
        expand(rules.download_manifest.output, exp=EXPERIMENTS)
