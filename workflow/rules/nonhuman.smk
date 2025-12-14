module NonHuman:
    snakefile:
        "nonhuman/Snakefile"
    config:
        {
            "output_dir": join(config["output_dir"], "nonhuman"),
            "logs_dir": join(config["logs_dir"], "nonhuman"),
            "benchmarks_dir": join(config["benchmarks_dir"], "nonhuman"),
            "samples": {
                "A_thaliana": {
                    # https://www.science.org/doi/full/10.1126/science.abi7489
                    "Col-CEN": {
                        "assembly": "https://github.com/schatzlab/Col-CEN/raw/refs/heads/main/v1.2/Col-CEN_v1.2.fasta.gz",
                        "hifi": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR621/003/ERR6210723/ERR6210723.fastq.gz",
                    },
                    # https://www.cell.com/molecular-plant/fulltext/S1674-2052%2822%2900181-2?dgcid=raven_jbs_etoc_email
                    "Col-PEK": {
                        "assembly": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/911/765/GCA_020911765.2_ASM2091176v2/GCA_020911765.2_ASM2091176v2_genomic.fna.gz",
                        "hifi": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18777598/SRR18777598",
                    },
                    # https://academic.oup.com/gpb/article/20/1/4/7230403
                    "Col-XJTU": {
                        "assembly": "https://download.cncb.ac.cn/gwh/Plants/Arabidopsis_thaliana_AT_GWHBDNP00000000.1/GWHBDNP00000000.1.genome.fasta.gz",
                        "hifi": "https://download.cncb.ac.cn/gsa/CRA004538/CRR302668/CRR302668.fastq.gz",
                    },
                }
            },
        }


use rule * from NonHuman as nonhuman_*
