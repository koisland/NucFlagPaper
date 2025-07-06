from os.path import join


FA = config["fa"]
SEED_OPTS = config["seed_opts"]
HAPS = ["MATERNAL", "PATERNAL"]

wildcard_constraints:
    seed="|".join(str(seed) for seed in SEED_OPTS.keys()),
    hap="|".join(HAPS),

rule compile_misasim:
    input:
        src="misasim/misasim",
    output:
        outbin=join(
            "results", "misasim", "release", "misasim"
        ),
    params:
        output_dir=lambda wc, output: os.path.dirname(os.path.dirname(str(output))),
    conda:
        "env.yaml"
    log:
        "logs/misasim/compile_misasim.log",
    shell:
        """
        log_file=$(realpath {log})
        output_dir=$(realpath {params.output_dir})
        cd {input}
        cargo build --release --target-dir ${{output_dir}} &> ${{log_file}}
        """

rule write_misasim_types:
    output:
        fa=join("results", "misasim", "simulated", "params.tsv"),
    run:
        with open(str(output), "wt") as fh:
            print("seed\tmtype\tnum\tmax_length", file=fh)
            for seed, (mtype, num, length) in SEED_OPTS.items():
                print(f"{seed}\t{mtype}\t{num}\t{length}", file=fh)


rule generate_misassemblies:
    input:
        bn=rules.compile_misasim.output,
        fa=FA,
    output:
        fa=join("results", "misasim", "simulated", "{seed}.fa"),
        bed=join("results", "misasim", "simulated", "{seed}.bed")
    params:
        mtype=lambda wc: SEED_OPTS[int(wc.seed)][0],
        num=lambda wc: SEED_OPTS[int(wc.seed)][1],
        length_arg=lambda wc: (
            f"-l {SEED_OPTS[int(wc.seed)][2]}"
            if SEED_OPTS[int(wc.seed)][0] != "break"
            else
            ""
        ),
        # Group by chromosome. Choose one hap to get misassembled.
        group_by=r"^(.*?)_.*?$"
    log:
        "logs/patch/generate_misassemblies_{seed}.log",
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
        join("results", "misasim", "simulated", "{seed}_{hap}.fa")
    conda:
        "env.yaml"
    shell:
        """
        samtools faidx {input.fa}
        seqtk subseq {input.fa} <(grep {wildcards.hap} "{input.fa}.fai" | cut -f 1) > {output}
        """

rule misasim_all:
    input:
        # Then generate misassemblies and check them with nucflag and flagger.
        expand(rules.generate_misassemblies.output, seed=SEED_OPTS.keys()),
        rules.write_misasim_types.output,
        expand(rules.split_asm_misasim.output, seed=SEED_OPTS.keys(), hap=HAPS),
    default_target: True
