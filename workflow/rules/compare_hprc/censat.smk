df_cens = pl.read_csv("data/compare_hprc/censat_hprc_r2_v1.0.index.csv")


rule download_r2_cens:
    output:
        join(OUTPUT_DIR, "data", "{sm}_R2_cens.bed"),
    params:
        uris=lambda wc: df_cens.filter(pl.col("sample_id").eq(pl.lit(wc.sm)))[
            "location"
        ].to_list(),
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        for uri in {params.uris}; do
            aws s3 cp ${{uri}} - | grep active_hor | awk -v OFS="\\t" '{{ print $1, $2, $3, "active_hor", $4}}' >> {output}
        done
        """


rule generate_cens_bed:
    input:
        fai=lambda wc: expand(
            rules.download_assemblies.output.asm_fai, sm=wc.sm, release="R2"
        ),
        censat=rules.download_r2_cens.output,
    output:
        join(OUTPUT_DIR, "data", "{sm}_R2_cens_other.bed"),
    conda:
        "../../envs/tools.yaml"
    shell:
        """
        cat \
            <(awk -v OFS="\\t" '{{ print $0, "cens"}}' {input.censat}) \
            <(
                bedtools subtract \
                    -a <(awk -v OFS="\\t" '{{ print $1, 0, $2 }}' {input.fai}) \
                    -b {input.censat} | \
                awk -v OFS="\\t" '{{ print $0, "other"}}'
            ) | \
        sort -k1,1 -k2,2n > {output}
        """


rule status_count_by_region:
    input:
        annot=rules.generate_cens_bed.output,
        calls=lambda wc: expand(
            rules.run_nucflag.output,
            sm=wc.sm,
            release="R2",
        ),
    output:
        join(OUTPUT_DIR, "censat", "{sm}_R2_status_{metric}_cens.bed"),
    conda:
        # This is a development branch I wrote specifically per Glennis's request.
        # Installed into this environment so not reproducible yet. Sorry.
        # https://github.com/logsdon-lab/NucFlag/commit/f3243363261366d03328cfc1bca7be060951d08b
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    shell:
        """
        nucflag status \
        -i {input.calls} \
        -b {input.annot} \
        -m {wildcards.metric} \
        -g name > {output}
        """


rule merge_statuses_and_label:
    input:
        expand(rules.status_count_by_region.output, sm=DATA.keys(), metric="{metric}"),
    output:
        join(OUTPUT_DIR, "censat", "all_R2_status_{metric}_cens.bed"),
    shell:
        """
        awk -v OFS="\\t" '{{
            if ($1 == "group") {{ next; }};
            match(FILENAME, "([^/]+)_R2", sms);
            print sms[1], $0
        }}' {input} > {output}
        """


rule censat_status_all:
    input:
        expand(rules.download_r2_cens.output, sm=DATA.keys()),
        expand(
            rules.status_count_by_region.output,
            sm=DATA.keys(),
            metric=["count", "length"],
        ),
        expand(rules.merge_statuses_and_label.output, metric=["count", "length"]),
