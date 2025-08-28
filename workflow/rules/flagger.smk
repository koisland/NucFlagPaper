import json
from os.path import join, dirname


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
            Literal["bias_annot_dir"]: str,
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


rule make_annot_json:
    input:
        wg_bed=rules.create_wg_bed.output.bed,
        annot_dir=lambda wc: SAMPLES[wc.sm].get("bias_annot_dir", []),
    output:
        annot_json=join(OUTPUT_DIR, "{sm}_annot.json"),
    run:
        import os
        import json
        import glob

        annot_beds = (
            glob.glob(os.path.join(input.annot_dir, "*.bed")) if input.annot_dir else []
        )
        annot_map = {
            os.path.splitext(os.path.basename(file))[0]: file for file in annot_beds
        }
        annot_map["whole_genome"] = input.wg_bed

        with open(output.annot_json, "wt") as fh:
            json.dump(annot_map, fh, indent=4)


rule bam_to_cov:
    input:
        bam=lambda wc: SAMPLES[wc.sm]["bam"],
        annot_json=rules.make_annot_json.output,
    output:
        cov=join(OUTPUT_DIR, "{sm}.cov"),
    params:
        baseline_annotation="whole_genome",
        run_bias_detection=lambda wc: (
            "--runBiasDetection" if SAMPLES[wc.sm].get("bias_annot_dir") else ""
        ),
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
        # Convert bam to cov.gz with bam2cov program
        bam2cov --bam {input.bam} \
        --output {output.cov} \
        --annotationJson {input.annot_json} \
        --threads {threads} \
        --baselineAnnotation {params.baseline_annotation} {params.run_bias_detection} &> {log}
        """


rule run_flagger:
    input:
        alpha=lambda wc: SAMPLES[wc.sm]["alpha"],
        cov=rules.bam_to_cov.output.cov,
    output:
        os.path.join(OUTPUT_DIR, "{sm}", "final_flagger_prediction.bed"),
    log:
        join(LOG_DIR, "run_flagger_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "run_flagger_{sm}.txt")
    resources:
        mem=MEM,
    singularity:
        "docker://mobinasri/flagger:v1.1.0"
    params:
        output_dir=lambda wc, output: dirname(output[0]),
    threads: THREADS
    shell:
        """
        mkdir -p {params.output_dir}
        hmm_flagger \
            --input {input.cov} \
            --outputDir {params.output_dir}  \
            --alphaTsv {input.alpha} \
            --labelNames Err,Dup,Hap,Col \
            --threads {threads} &> {log}
        """


rule flagger:
    input:
        expand(rules.run_flagger.output, sm=SAMPLES.keys()),
    default_target: True
