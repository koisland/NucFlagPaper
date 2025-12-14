

rule align_element:
    input:
        ref=rules.download_curated_asm.output.fa,
        reads=rules.download_element.output,
    output:
        bam=join(OUTPUT_DIR, "HG002_{version}_element.bam"),
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
        samtools fixmate -m - | \
        samtools sort -O bam -o {output.bam}_tmp ;}} 2>> {log}

        # mark duplicates and remove them
        # collect primary alignments
        {{ samtools markdup -r -@{threads} {output.bam}_tmp | \
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
        bed=join(OUTPUT_DIR, "HG002_{version}_lcr.bed"),
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
        bed=rules.find_all_lcr_regions.output,
        calls=expand(
            rules.versioned_check_asm_nucflag.output.misassemblies,
            sm="HG002_{version}",
        ),
    log:
        join("logs", "curated", "label_lcr_regions_and_pos_HG002_{version}.log"),
    shell:
        """
        """


# rule check_element_bam:
#     input:
#         ref=expand(rules.download_curated_asm.output.fa, version=f"v1.0.1")[0],
#         bam="",
#         bai="",
#     output:
#         misassemblies="",
#     threads: 12
#     shell:
#         """
#         nucflag call \
#         -i {input.bam} \
#         -o {output.misassemblies} \
#         -t {params.threads} \
#         -p {params.threads}
#         """


rule validate_homopolymers_all:
    input:
        expand(rules.align_element.output, version=ASSEMBLIES.keys()),
        expand(rules.find_all_lcr_regions.output, version=ASSEMBLIES.keys()),
    default_target: True
