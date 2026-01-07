from os.path import join

"""
{
    Literal["output_dir"]: str,
    Literal["logs_dir"]: str,
    Literal["benchmarks_dir"]: str,
    Literal["threads"]: int,
    Literal["mem"]: str,
    Literal["samples"]: [
        {
            Literal["name"]: str,
            Literal["asm_fa"]: str,
            Literal["bam"]: str,
        }
    ]
}
"""
SAMPLES = {sm["name"]: sm for sm in config["samples"]}
OUTPUT_DIR = config["output_dir"]
BENCHMARK_DIR = config["benchmarks_dir"]
LOG_DIR = config["logs_dir"]
THREADS = config.get("threads", 12)
MEM = config.get("mem", "250GB")


rule run_deepvariant:
    input:
        ref=lambda wc: SAMPLES[wc.sm]["asm_fa"],
        bam=lambda wc: SAMPLES[wc.sm]["bam"],
    output:
        vcf=join(OUTPUT_DIR, "{sm}.vcf.gz"),
    singularity:
        "docker://google/deepvariant:1.9.0"
    params:
        model="PACBIO",
    log:
        join(LOG_DIR, "run_deepvariant_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "run_deepvariant_{sm}.txt")
    resources:
        mem=MEM,
    threads: THREADS
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="{params.model}" \
        --vcf_stats_report=true \
        --ref="{input.ref}" \
        --reads="{input.bam}" \
        --output_vcf="{output.vcf}" \
        --num_shards={threads} 2> {log}
        """


"""
Filter VCF calls to those that PASS
Also consider GQ?
* https://github.com/google/deepvariant/issues/503
* https://github.com/google/deepvariant/issues/278

DeepPolisher defines a filter threshold based on variant type.
* https://pmc.ncbi.nlm.nih.gov/articles/PMC12212083/#SM1
* From Q100 manuscript:
    * 'GQ > 20 for 1bp insertions, GQ > 12 for 1bp deletions, and GQ > 5 for all other edit sizes, as recommended in Mastoras et al.48'
"""


rule filter_vcf_calls:
    input:
        vcf=rules.run_deepvariant.output.vcf,
    output:
        bed=join(OUTPUT_DIR, "{sm}.bed"),
    shell:
        """
        zcat -f {input.vcf} | grep -v "#" | awk -v OFS="\\t" '{{
            # TODO: Filter by GQ?
            if ($7 != "PASS") {{ next }};
            print $1, $2, $2 + length($4), $4"-"$5, 0, $3, $2, $2 + length($4), "0,0,0"
        }}' > {output}
        """


rule deepvariant:
    input:
        expand(rules.run_deepvariant.output, sm=SAMPLES.keys()),
        expand(rules.filter_vcf_calls.output, sm=SAMPLES.keys()),
    default_target: True
