

rule generate_index:
    input:
        ref=rules.download_curated_asm.output.fa,
    output:
        multiext(
            rules.download_curated_asm.output.fa,
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    conda:
        "../../envs/curated.yaml"
    log:
        join("logs", "curated", "index_HG002_{version}.log"),
    shell:
        """
        bwa index {input.ref} &> {log}
        """


rule align_element:
    input:
        ref=rules.download_curated_asm.output.fa,
        idx=rules.generate_index.output,
        reads=rules.download_element.output,
    output:
        bam=temp(join(OUTPUT_DIR, "HG002_{version}_element.tmp.bam")),
    threads: 12 + 1 + 1
    conda:
        "../../envs/curated.yaml"
    log:
        join("logs", "curated", "align_HG002_{version}_element.log"),
    shell:
        """
        # From Q100 manuscript (Arang's T2T-polish bwa pipeline workflow)
        # https://pmc.ncbi.nlm.nih.gov/articles/instance/12458380/bin/media-1.pdf
        {{ bwa mem -t {threads} {input.ref} {input.reads} | \
        samtools fixmate -m - - | \
        samtools sort -O bam -o {output.bam} ;}} 2>> {log}
        """


rule alignment_element_remove_dup:
    input:
        bam=rules.align_element.output,
    output:
        bam=join(OUTPUT_DIR, "HG002_{version}_element.bam"),
        bai=join(OUTPUT_DIR, "HG002_{version}_element.bam.bai"),
    threads: 12
    conda:
        "../../envs/curated.yaml"
    log:
        join("logs", "curated", "remove_dup_HG002_{version}_element.log"),
    shell:
        """
        # mark duplicates and remove them
        # collect primary alignments
        {{ samtools markdup -r -@{threads} {input.bam} - | \
        samtools view -F0x100 -hb -o {output.bam} ;}} 2>> {log}
        samtools index -@{threads} {output.bam} 2>> {log}
        """


"""
Where are all homopolymers/simple_repeats?
"""


rule find_all_lcr_regions:
    input:
        ref=rules.download_curated_asm.output.fa,
    output:
        bed=join(OUTPUT_DIR, "homopolymers", "HG002_{version}_lcr.bed"),
    conda:
        "../../envs/curated.yaml"
    log:
        join("logs", "curated", "find_all_lcr_regions_HG002_{version}.log"),
    shell:
        """
        longdust {input.ref} > {output.bed} 2> {log}
        """


"""
What distribution of nucflag calls have abnormal signals
Where do these errors fall in? Based on observation that always at start of hp run
"""


rule label_lcr_regions_and_pos:
    input:
        script="workflow/scripts/curated/label_homopolymers.py",
        ref=rules.download_curated_asm.output.fa,
        bed=rules.find_all_lcr_regions.output,
        calls=expand(
            rules.versioned_check_asm_nucflag.output.misassemblies,
            sm="HG002_{version}",
        ),
    output:
        join(OUTPUT_DIR, "homopolymers", "HG002_{version}_homopolymers.bed"),
    conda:
        "../../envs/curated.yaml"
    log:
        join("logs", "curated", "label_lcr_regions_and_pos_HG002_{version}.log"),
    shell:
        """
        # bgzip weird
        zcat {input.ref} > {input.ref}.tmp
        samtools faidx {input.ref}.tmp
        python {input.script} -i {input.bed} -r {input.ref}.tmp -c {input.calls}
        rm -f {input.ref}.tmp {input.ref}.tmp.fai
        """


rule versioned_element_check_asm_nucflag:
    input:
        bam=rules.alignment_element_remove_dup.output.bam,
        fa=rules.download_curated_asm.output.fa,
        fai=rules.download_curated_asm.output.fa + ".fai",
        config=CONFIG_ELEMENT,
    output:
        misassemblies=os.path.join(
            OUTPUT_DIR,
            "HG002_{version}_element_calls.bed",
        ),
        pileup_dir=directory(join(OUTPUT_DIR, "HG002_{version}_element_pileup")),
    params:
        config=CONFIG_ELEMENT,
    conda:
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    threads: 12
    resources:
        threads=12,
        mem="50GB",
    log:
        os.path.join(LOGS_DIR, "run_nucflag_element_HG002_{version}.log"),
    benchmark:
        os.path.join(BMKS_DIR, "run_nucflag_element_HG002_{version}.tsv")
    shell:
        """
        nucflag call \
        -i {input.bam} \
        -f {input.fa} \
        -c {input.config} \
        -o {output.misassemblies} \
        -t {resources.threads} \
        --output_pileup_dir {output.pileup_dir} \
        --add_pileup_data mapq \
        -p {threads} 2> {log}
        """


"""
Get consensus of callset. Element data is less vulnerable to homopolymer expansions and contractions.
We can validate the hifi calls this way.
"""


rule get_consensus_nucflag:
    input:
        calls_element=rules.versioned_element_check_asm_nucflag.output,
        calls_hifi=expand(
            rules.versioned_check_asm_nucflag.output.misassemblies,
            sm="HG002_{version}",
        ),
    output:
        consensus=os.path.join(
            OUTPUT_DIR,
            "HG002_{version}_hifi-element_consensus_calls.bed",
        ),
    params:
        filter_call_pattern="homopolymer|dinucleotide|simple_repeat",
    conda:
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    log:
        os.path.join(LOGS_DIR, "run_nucflag_element_HG002_{version}.log"),
    shell:
        """
        nucflag consensus \
        -i <(grep -P "{params.filter_call_pattern}" {input.calls_hifi}) \
        <(grep -P "{params.filter_call_pattern}" {input.calls_element}) > {output} 2> {log}
        """


# Merge MAPQ bigWigs and convert to bedGraph.
rule merge_element_mapq_bigwigs:
    input:
        bw_dir_element=rules.versioned_element_check_asm_nucflag.output.pileup_dir,
    output:
        bw=os.path.join(
            OUTPUT_DIR,
            "HG002_{version}_element_mapq.bw",
        ),
        bg=os.path.join(
            OUTPUT_DIR,
            "HG002_{version}_element_mapq.bg.gz",
        ),
    params:
        bg_interm=lambda wc, output: output.bg.replace(".gz", ""),
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        bigwigmerge -l <(find {input.bw_dir_element}/*.bw) {output.bw}
        bigwigtobedgraph {output.bw} {params.bg_interm} && bgzip {params.bg_interm}
        """


# TODO: Summarize/plot checking against low MAPQ regions. Histogram of valid hits


rule validate_homopolymers_all:
    input:
        expand(rules.alignment_element_remove_dup.output, version=ASSEMBLIES.keys()),
        expand(
            rules.versioned_element_check_asm_nucflag.output, version=ASSEMBLIES.keys()
        ),
        expand(rules.find_all_lcr_regions.output, version=ASSEMBLIES.keys()),
        expand(rules.get_consensus_nucflag.output, version=ASSEMBLIES.keys()),
        expand(rules.merge_element_mapq_bigwigs.output, version=ASSEMBLIES.keys()),
    default_target: True
