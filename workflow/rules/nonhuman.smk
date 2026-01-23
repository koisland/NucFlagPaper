module NonHuman:
    snakefile:
        "nonhuman/Snakefile"
    config:
        {
            "output_dir": join(config["output_dir"], "nonhuman"),
            "logs_dir": join(config["logs_dir"], "nonhuman"),
            "benchmarks_dir": join(config["benchmarks_dir"], "nonhuman"),
            "reference": "Col-PEK",
            "samples": {
                # https://www.science.org/doi/full/10.1126/science.abi7489
                "Col-CEN": {
                    "assembly": "https://github.com/schatzlab/Col-CEN/raw/refs/heads/main/v1.2/Col-CEN_v1.2.fasta.gz",
                    "hifi": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR621/003/ERR6210723/ERR6210723.fastq.gz",
                    # Two replicates are available but most reads are duplicate.
                    "ont": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR553/005/ERR5530735/ERR5530735.fastq.gz",
                    "config": config["config"]["hifi_curated"],
                    "bed_centromeres": config["nonhuman"]["bed_centromeres"]["Col-CEN"],
                },
                # https://www.cell.com/molecular-plant/fulltext/S1674-2052%2822%2900181-2?dgcid=raven_jbs_etoc_email
                "Col-PEK": {
                    "assembly": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/911/765/GCA_020911765.2_ASM2091176v2/GCA_020911765.2_ASM2091176v2_genomic.fna.gz",
                    "hifi": "https://download.cncb.ac.cn/gsa/CRA005381/CRR339253/CRR339253.fastq.gz",
                    "ont": "https://download.cncb.ac.cn/gsa/CRA005345/CRR338453/CRR338453.fastq.gz",
                    "config": config["config"]["hifi_curated"],
                    "bed_centromeres": config["nonhuman"]["bed_centromeres"]["Col-PEK"],
                },
                # https://academic.oup.com/gpb/article/20/1/4/7230403
                "Col-XJTU": {
                    "assembly": "https://download.cncb.ac.cn/gwh/Plants/Arabidopsis_thaliana_AT_GWHBDNP00000000.1/GWHBDNP00000000.1.genome.fasta.gz",
                    "hifi": "https://download.cncb.ac.cn/gsa/CRA004538/CRR302668/CRR302668.fastq.gz",
                    "ont": "https://download.cncb.ac.cn/gsa/CRA004538/CRR302667/CRR302667.fastq.gz",
                    "config": config["config"]["hifi_curated"],
                    "bed_centromeres": config["nonhuman"]["bed_centromeres"][
                        "Col-XJTU"
                    ],
                },
            },
        }


use rule * from NonHuman as nonhuman_*
