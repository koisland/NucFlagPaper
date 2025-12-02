
rule download_element:
    output:
        reads_1=join(OUTPUT_DIR, "data", "reads", "ASHG-C063_R1.fastq.gz"),
        reads_2=join(OUTPUT_DIR, "data", "reads", "ASHG-C063_R2.fastq.gz"),
    params:
        uri_1="s3://human-pangenomics/T2T/scratch/HG002/sequencing/element/trio/HG002/ins1000/ASHG-C063_R1.fastq.gz",
        uri_2="s3://human-pangenomics/T2T/scratch/HG002/sequencing/element/trio/HG002/ins1000/ASHG-C063_R2.fastq.gz",
    shell:
        """
        aws s3 cp {params.uri_1} {output.reads_1}
        aws s3 cp {params.uri_2} {output.reads_2}
        """


rule align_element:
    input:
        ref=expand(rules.download_curated_asm.output.fa, version=f"v1.0.1")[0],
        reads="",
    output:
        bam="",
    threads: 12 + 1 + 1
    shell:
        """
        # From Q100 manuscript (Arang's T2T-polish bwa pipeline workflow)
        # https://pmc.ncbi.nlm.nih.gov/articles/instance/12458380/bin/media-1.pdf
        bwa mem -t {threads} {input.ref} {input.reads} | \
        samtools fixmate -m - | \
        samtools sort -O bam -o {output.bam}_tmp
        # samtools index {threads} {output.bam}_tmp

        # mark duplicates and remove them
        # collect primary alignments
        samtools markdup -r -@{threads} {output.bam}_tmp | \
        samtools view -F0x100 -hb -o {output.bam}
        samtools index -@{threads} {output.bam}
        """


rule check_bam:
    input:
        ref=expand(rules.download_curated_asm.output.fa, version=f"v1.0.1")[0],
        bam="",
        bai="",
    output:
        misassemblies="",
    threads: 12
    shell:
        """
        nucflag call \
        -i {input.bam} \
        -o {output.misassemblies} \
        -t {params.threads} \
        -p {params.threads}
        """
