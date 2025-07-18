from os.path import join

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
MEM = config.get("mem", "50GB")


rule run_inspector:
    input:
        bam=lambda wc: SAMPLES[wc.sm]["bam"],
        fa=lambda wc: SAMPLES[wc.sm]["asm_fa"],
    output:
        output_dir=directory(join(OUTPUT_DIR, "{sm}")),
    params:
        datatype=lambda wc: SAMPLES[wc.sm]["datatype"],
    log:
        join(LOG_DIR, "run_inspector_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "run_inspector_{sm}.txt"),
    resources:
        mem=MEM,
    threads:
        THREADS
    conda:
        "../../envs/misasim.yaml"
    shell:
        """
        # https://github.com/Maggi-Chen/Inspector/issues/10
        # Symlink BAMs
        ln -s $(realpath {input.bam}) {output.output_dir}/read_to_contig.bam
        ln -s $(realpath {input.bam}.bai) {output.output_dir}/read_to_contig.bam.bai

        # Skip read alignment
        inspector.py \
        -t {threads} \
        -c {input.fa} \
        -r {input.fa} \
        -o {output.output_dir} \
        --datatype {params.datatype} \
        --skip_read_mapping &> {log}
        """


rule inspector:
    input:
        expand(rules.run_inspector.output, sm=SAMPLES.keys()),
    default_target: True
