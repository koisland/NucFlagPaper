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
            Literal["alpha"]: str,
        }
    ]
}
"""
SAMPLES = {sm["name"]: sm for sm in config["samples"]}
OUTPUT_DIR = config["output_dir"]
BENCHMARK_DIR = config["benchmarks_dir"]
LOG_DIR = config["logs_dir"]
THREADS = config.get("threads", 12)
MEM = config.get("mem", "100GB")


wildcard_constraints:
    sm="|".join(SAMPLES.keys()),


rule create_wg_bed:
    input:
        fai=lambda wc: SAMPLES[wc.sm]["asm_fa"] + ".fai",
    output:
        bed=join(OUTPUT_DIR, "{sm}_wg.bed"),
    shell:
        """
        awk -v OFS="\\t" '{{ print $1, 0, $2 }}' {input.fai} > {output.bed}
        """


rule bam_to_cov:
    input:
        wg_bed=rules.create_wg_bed.output,
        bam=lambda wc: SAMPLES[wc.sm]["bam"],
    output:
        annot_json=join(OUTPUT_DIR, "{sm}_annot.json"),
        cov=join(OUTPUT_DIR, "{sm}.cov"),
    params:
        baseline_annotation="whole_genome",
    log:
        join(LOG_DIR, "bam_to_cov_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "bam_to_cov_{sm}.txt")
    singularity:
        "docker://mobinasri/flagger:v1.1.0"
    resources:
        mem=MEM,
    threads: THREADS
    shell:
        """
        # Put the path to the whole-genome bed file in a json file
        echo "{{" > {output.annot_json}
        echo '"whole_genome": "{input.wg_bed}"' >> {output.annot_json}
        echo "}}" >> {output.annot_json}

        # Convert bam to cov.gz with bam2cov program
        bam2cov --bam {input.bam} \
        --output {output.cov} \
        --annotationJson {output.annot_json} \
        --threads {threads} \
        --baselineAnnotation {params.baseline_annotation} &> {log}
        """


rule run_flagger:
    input:
        alpha=lambda wc: SAMPLES[wc.sm]["alpha"],
        cov=rules.bam_to_cov.output.cov,
    output:
        output_dir=directory(os.path.join(OUTPUT_DIR, "{sm}")),
    log:
        join(LOG_DIR, "run_flagger_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "run_flagger_{sm}.txt")
    resources:
        mem=MEM,
    singularity:
        "docker://mobinasri/flagger:v1.1.0"
    threads: THREADS
    shell:
        """
        mkdir -p {output.output_dir}
        hmm_flagger \
            --input {input.cov} \
            --outputDir {output.output_dir}  \
            --alphaTsv {input.alpha} \
            --labelNames Err,Dup,Hap,Col \
            --threads {threads} &> {log}
        """


rule flagger:
    input:
        expand(rules.run_flagger.output, sm=SAMPLES.keys()),
    default_target: True
