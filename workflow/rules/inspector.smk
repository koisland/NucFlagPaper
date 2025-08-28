from os.path import join, abspath

"""
{
    Literal["output_dir"]: str,
    Literal["logs_dir"]: str,
    Literal["benchmarks_dir"]: str,
    Literal["samples"]: [
        {
            Literal["name"]: str,
            Literal["asm_fa"]: str,
            Literal["datatype"]: str,
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


rule run_inspector:
    input:
        bam=lambda wc: SAMPLES[wc.sm]["bam"],
        fa=lambda wc: SAMPLES[wc.sm]["asm_fa"],
    output:
        output_dir=directory(abspath(join(OUTPUT_DIR, "{sm}"))),
    params:
        datatype=lambda wc: SAMPLES[wc.sm]["datatype"],
        # Keep in params as if any core dump will trigger rerun.
        inspector_dir="workflow/scripts/Inspector",
    log:
        abspath(join(LOG_DIR, "run_inspector_{sm}.log")),
    benchmark:
        join(BENCHMARK_DIR, "run_inspector_{sm}.txt")
    resources:
        mem=MEM,
    threads: THREADS
    conda:
        "../envs/inspector.yaml"
    shell:
        """
        # https://github.com/Maggi-Chen/Inspector/issues/10
        # Symlink BAMs
        mkdir -p {output.output_dir}
        bam=$(realpath {input.bam})
        fa=$(realpath {input.fa})
        ln -s ${{bam}} {output.output_dir}/read_to_contig.bam
        ln -s ${{bam}}.bai {output.output_dir}/read_to_contig.bam.bai

        # Skip read alignment
        cd {params.inspector_dir}
        python inspector.py \
        -t {threads} \
        -c ${{fa}} \
        -r ${{fa}} \
        -o {output.output_dir} \
        --datatype {params.datatype} \
        --skip_read_mapping &> {log}
        """


# All calls are significant already.
rule merge_calls:
    input:
        inspector_dir=rules.run_inspector.output,
    output:
        bed=join(OUTPUT_DIR, "{sm}.bed"),
    shell:
        """
        sort -k1,1 -k2,2n \
        <(awk -v OFS="\t" 'NR > 1 {{ print $1, $2, $3, $5}}' {input.inspector_dir}/structural_error.bed) \
        <(awk -v OFS="\t" 'NR > 1 {{ print $1, $2, $3, $8}}' {input.inspector_dir}/small_scale_error.bed) > {output}
        """


rule inspector:
    input:
        expand(rules.run_inspector.output, sm=SAMPLES.keys()),
        expand(rules.merge_calls.output, sm=SAMPLES.keys()),
    default_target: True
