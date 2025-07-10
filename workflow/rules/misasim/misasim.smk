from os.path import join


HAPS = ["MATERNAL", "PATERNAL"]

# {str: {Literal["fa"]: str, Literal["seed_opts"]: dict[int, tuple[str, int, int]]}}
SAMPLE_OPTS = config["samples"]
OUTPUT_DIR = config.get("output_dir", "results/misasim")
LOG_DIR = config.get("log_dir", "logs/misasim")
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
        output_dir=lambda wc, output: os.path.dirname(os.path.dirname(str(output))),
    log:
        join(LOG_DIR, "compile_misasim.log"),
    # TODO: Need cargo.
    shell:
        """
        log_file=$(realpath {log})
        output_dir=$(realpath {params.output_dir})
        cd {input}
        cargo build --release --target-dir ${{output_dir}} &> ${{log_file}}
        """


rule write_misasim_types:
    output:
        fa=join(OUTPUT_DIR, "{sm}", "params.tsv"),
    run:
        with open(str(output), "wt") as fh:
            print("seed\tmtype\tnum\tmax_length", file=fh)
            for seed, (mtype, num, length) in SAMPLE_OPTS[wildcards.sm][
                "seed_opts"
            ].items():
                print(f"{seed}\t{mtype}\t{num}\t{length}", file=fh)


rule generate_misassemblies:
    input:
        bn=rules.compile_misasim.output,
        fa=lambda wc: SAMPLE_OPTS[wc.sm]["fa"],
    output:
        fa=join(OUTPUT_DIR, "{sm}", "{seed}.fa"),
        bed=join(OUTPUT_DIR, "{sm}", "{seed}.bed"),
    params:
        mtype=lambda wc: SAMPLE_OPTS[wc.sm]["seed_opts"][int(wc.seed)][0],
        num=lambda wc: SAMPLE_OPTS[wc.sm]["seed_opts"][int(wc.seed)][1],
        length_arg=lambda wc: (
            f"-l {SAMPLE_OPTS[wc.sm]['seed_opts'][int(wc.seed)][2]}"
            if SAMPLE_OPTS[wc.sm]["seed_opts"][int(wc.seed)][0] != "break"
            else ""
        ),
        # Group by chromosome. Choose one hap to get misassembled.
        group_by=r"^(.*?)_.*?$",
    log:
        join(LOG_DIR, "generate_misassemblies_{sm}_{seed}.log"),
    shell:
        """
        {input.bn} {params.mtype} \
        -n {params.num} \
        {params.length_arg} \
        -s {wildcards.seed} \
        -i {input.fa} \
        -o {output.fa} \
        -g "{params.group_by}" \
        -b {output.bed} 2> {log}
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
        expand(rules.write_misasim_types.output, sm=SAMPLES),
        expand(rules.split_asm_misasim.output, zip, sm=SAMPLES, seed=SEEDS, hap=HAPS),
    default_target: True
