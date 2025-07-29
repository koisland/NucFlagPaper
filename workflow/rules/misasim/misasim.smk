import json
from os.path import join, dirname


HAPS = ["MATERNAL", "PATERNAL"]

# {str: {Literal["fa"]: str, Literal["seed_opts"]: dict[int, dict[str, Any]]}}
SAMPLE_OPTS = config["samples"]
OUTPUT_DIR = config.get("output_dir", "results/misasim")
LOG_DIR = config.get("log_dir", "logs/misasim")
GROUP_BY = config.get("group_by", r"^(.*?)_(.*?)$")
SAMPLES = []
SEEDS = []
for sm, opts in SAMPLE_OPTS.items():
    for seed, seed_opts in opts["seed_opts"].items():
        SAMPLES.append(sm)
        SEEDS.append(str(seed))


wildcard_constraints:
    sm="|".join(SAMPLES),
    seed="|".join(SEEDS),
    hap="|".join(HAPS),


rule compile_misasim:
    input:
        src="workflow/rules/misasim/misasim",
    output:
        outbin=join(OUTPUT_DIR, "release", "misasim"),
    params:
        output_dir=lambda wc, output: dirname(dirname(str(output))),
    log:
        join(LOG_DIR, "compile_misasim.log"),
    shell:
        """
        log_file=$(realpath {log})
        output_dir=$(realpath {params.output_dir})
        cd {input}
        cargo build --release --target-dir ${{output_dir}} &> ${{log_file}}
        """


rule generate_misassemblies:
    input:
        bn=rules.compile_misasim.output,
        fa=lambda wc: SAMPLE_OPTS[wc.sm]["fa"],
    output:
        fa=join(OUTPUT_DIR, "{sm}", "{seed}.fa"),
        cfg=join(OUTPUT_DIR, "{sm}", "{seed}.json"),
        bed=join(OUTPUT_DIR, "{sm}", "{seed}.bed"),
    params:
        # Group by chromosome. Choose one hap to get misassembled.
        # r"^(.*?)_.*?$"
        # No group. All haps get misassembled.
        # r"^(.*?)_(.*?)$"
        group_by=f'-g "{GROUP_BY}"' if GROUP_BY else "",
        # JSON str with seed options.
        # Must be array. ([{"mtype": ..., "number": ..., "length": ...}])
        cfg=lambda wc: json.dumps([SAMPLE_OPTS[wc.sm]["seed_opts"][int(wc.seed)]]),
    log:
        join(LOG_DIR, "generate_misassemblies_{sm}_{seed}.log"),
    shell:
        """
        echo '{params.cfg}' > {output.cfg}
        {input.bn} multiple \
        -p {output.cfg} \
        -s {wildcards.seed} \
        -i {input.fa} \
        -o {output.fa} \
        {params.group_by} \
        -b {output.bed} 2> {log}
        """


rule index_misassemblies_fa:
    input:
        fa=rules.generate_misassemblies.output.fa,
    output:
        idx=join(OUTPUT_DIR, "{sm}", "{seed}.fa.fai"),
    conda:
        "../../envs/misasim.yaml"
    shell:
        """
        samtools faidx {input}
        """


rule split_asm_misasim:
    input:
        fa=rules.generate_misassemblies.output.fa,
    output:
        join(OUTPUT_DIR, "{sm}", "{seed}_{hap}.fa"),
    conda:
        "../../envs/misasim.yaml"
    shell:
        """
        samtools faidx {input.fa}
        seqtk subseq {input.fa} <(grep {wildcards.hap} "{input.fa}.fai" | cut -f 1) > {output}
        """


rule all:
    input:
        # Then generate misassemblies and check them with nucflag and flagger.
        expand(rules.generate_misassemblies.output, zip, sm=SAMPLES, seed=SEEDS),
        expand(rules.index_misassemblies_fa.output, zip, sm=SAMPLES, seed=SEEDS),
        expand(rules.split_asm_misasim.output, zip, sm=SAMPLES, seed=SEEDS, hap=HAPS),
    default_target: True
