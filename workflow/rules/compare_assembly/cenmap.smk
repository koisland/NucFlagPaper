
rule run_cenmap:
    input:
        asm=lambda wc: get_assembly(wc.sm),
    output:
        bed=join(
            OUTPUT_DIR,
            "cenmap",
            "{sm}",
            "{sm}",
            "results",
            "final",
            "bed",
            "all_AS-HOR_lengths.bed",
        ),
        fai=join(
            OUTPUT_DIR,
            "cenmap",
            "{sm}",
            "{sm}",
            "results",
            "2-concat_asm",
            "{sm}-asm-renamed-reort.fa.fai",
        ),
    params:
        output_dir=join(OUTPUT_DIR, "cenmap", "{sm}"),
        profile="workflow/profiles/lpc_all",
    conda:
        "../../envs/cenmap.yaml"
    log:
        join(LOG_DIR, "cenmap", "{sm}.log"),
    threads: 1
    shell:
        """
        cenmap \
        -i {input.asm} \
        -o {params.output_dir} \
        -j 20 \
        --workflow-profile {params.profile} \
        -s {wildcards.sm} 2> {log}
        """


rule create_hor_array_bed:
    input:
        bed=rules.run_cenmap.output.bed,
        fai=rules.run_cenmap.output.fai,
    output:
        join(
            OUTPUT_DIR,
            "cenmap",
            "{sm}",
            "{sm}",
            "results",
            "final",
            "bed",
            "all_AS-HOR_lengths_original_coords.bed",
        ),
    shell:
        """
        # Remove coordinates
        # Sort
        # Join with contig lengths
        # Reverse coordinates if needed and strip chrom and sample from name to get original contig name.
        awk -v OFS="\\t" '{{
            match($1, "^(.+):", arr);
            print arr[1], $2, $3
        }}' {input.bed} | \
        sort -k1,1 | \
        join - <(sort -k1,1 {input.fai} | cut -f 1,2) | \
        awk -v OFS="\\t" '{{
            if (match($1, "rc-")) {{ st=$4-$3; end=$4-$2 }} else {{ st=$2; end=$3 }};
            sub("{wildcards.sm}_", "", $1);
            sub("rc-", "", $1);
            sub("chr[0-9XY]+_", "", $1);
            print $1, st, end, $1
        }}' > {output}
        """


rule generate_status_bed:
    input:
        bed=rules.create_hor_array_bed.output,
        calls=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm="{sm}_hifi",
        ),
    output:
        join(
            OUTPUT_DIR,
            "cenmap",
            "{sm}",
            "{sm}",
            "results",
            "final",
            "bed",
            "nucflag_status.bed",
        ),
    conda:
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    shell:
        """
        nucflag status -i {input.calls} -b {input.bed} -g name > {output}
        """


rule cenmap_all:
    input:
        expand(rules.run_cenmap.output, sm=config["samples"].keys()),
        expand(rules.generate_status_bed.output, sm=config["samples"].keys()),
