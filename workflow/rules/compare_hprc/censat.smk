df_censat = pl.read_csv(
    "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/censat/censat_hprc_r2_v1.0.index.csv"
)


rule download_r2_censat:
    output:
        join(OUTPUT_DIR, "data", "{sm}_R2_censat.bed"),
    params:
        uris=lambda wc: df_censat.filter(pl.col("sample_id").eq(pl.lit(wc.sm)))[
            "location"
        ].to_list(),
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        for uri in {params.uris}; do
            aws s3 cp ${{uri}} - >> {output}
        done
        """


rule status_by_region:
    input:
        censat=rules.download_r2_censat.output,
        calls=lambda wc: expand(
            rules.run_nucflag.output,
            sm=wc.sm,
            release="R2",
        ),
    output:
        join(OUTPUT_DIR, "censat", "{sm}_R2_status_censat.bed"),
    conda:
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    shell:
        """
        nucflag status \
        -i {input.calls} \
        -b <(grep -v track {input.censat} | \
            awk -v OFS="\\t" '{{
                match($4, "^(.+)\\\\(", res);
                print $1, $2, $3, res[1] == "" ? $4 : res[1]
            }}') \
        -g name > {output}
        """


rule censat_status_all:
    input:
        expand(rules.download_r2_censat.output, sm=DATA.keys()),
        expand(rules.status_by_region.output, sm=DATA.keys()),
