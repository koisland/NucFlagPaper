ANNOT_COLORS = {
    "censat": "config/censat_colors.tsv",
    "segdups": "config/segdups_colors.tsv",
}
ANNOT_TITLE = {"censat": "Centromere satellites", "segdups": "Segmental duplications"}


wildcard_constraints:
    annot="|".join(ANNOT_COLORS.keys())


rule plot_compare_misassemblies:
    input:
        calls=[
            expand(
                rules.nf_denovo_check_asm_nucflag.output.misassemblies,
                sm=f"{sm}_hifi",
            )
            for sm in ASM_COLORS.keys()
        ],
    output:
        plot=multiext(
            join(OUTPUT_DIR, "nucflag", "all_calls_compared"), ".png", ".pdf", ".tsv"
        ),
    params:
        script="workflow/scripts/compare_assembly/cmp_all_calls.py",
        calls=lambda wc, input: " ".join(
            [
                f"<(grep -v chrY {call})" if "CHM13v2.0" in call else call
                for call in input.calls
            ]
        ),
        labels=" ".join(ASM_COLORS.keys()),
        colors=" ".join(ASM_COLORS.values()),
        output_prefix=lambda wc, output: splitext(output.plot[0])[0],
    conda:
        "../../envs/curated.yaml"
    shell:
        """
        python {params.script} -i {params.calls} -l {params.labels} -c {params.colors} -o {params.output_prefix}
        """


rule aggregate_call_num:
    input:
        calls=[
            expand(
                rules.nf_denovo_check_asm_nucflag.output.misassemblies,
                sm=f"{sm}_hifi",
            )
            for sm in ASM_COLORS.keys()
        ],
        lengths=rules.plot_compare_misassemblies.output[2],
    output:
        status=join(OUTPUT_DIR, "nucflag", "agg_errors.tsv"),
    shell:
        """
        awk -v OFS="\\t" 'FNR > 1 {{
            if ($4 == "#chrom") {{ next; }}
            match(FILENAME, "([^\\/]+)_hifi", arr);
            print $4, arr[1]
        }}' {input.calls} | \
        sort -k2,2 -k1,1 | \
        uniq -c | \
        awk -v OFS="\\t" '{{ print $2, $3, $1 }}' > {output.status}.tmp
        # Join to add based on type and assembly
        join -a 1 \
            <(sort -k1,1 -k2,2 {output.status}.tmp | awk -v OFS="\\t" '{{ print $1"~"$2, $3 }}') \
            <(awk 'NR > 1 {{ print $2"~"$1, $3 }}' {input.lengths} | sort -k1,1 -k2,2) | \
        awk -v OFS="\\t" '{{
            split($1, arr, "~");
            $1=$1;
            print arr[1], arr[2], $2, $3
        }}' > {output.status}
        rm -f {output.status}.tmp
        """


def get_annotations(wc):
    if wc.sm == "CHM13v2.0":
        return getattr(rules.download_annotations.output, wc.annot)
    else:
        return expand(
            getattr(rules.liftover_annotations.output, wc.annot),
            ref="CHM13v2.0",
            sm=wc.sm,
        )


def get_awk_cmd_by_annot(wc) -> str:
    if wc.annot == "segdups":
        return """
            awk -v OFS="\\t" '{{
                sim=$24;
                if (sim < 0.9) {{
                    lbl = "Less than 90%"
                }} else if (sim >= 0.9 && sim < 0.98) {{
                    lbl="90\\% - 98\\%"
                }} else if (sim >= 0.98 && sim < 0.99 ) {{
                    lbl="98\\% - 99\\%"
                }} else {{
                    lbl="Greater than 99%"
                }};
                print $1, $2, $3, lbl
            }}'"""
    elif wc.annot == "censat":
        return (
            """awk -v OFS="\\t" '{{ split($4, arr, "_"); print $1, $2, $3, arr[1]}}'"""
        )
    else:
        raise ValueError(f"Invalid annotation {wc.annot}")


rule generate_status_by_annot:
    input:
        annot=get_annotations,
        calls=expand(
            rules.nf_denovo_check_asm_nucflag.output.misassemblies,
            sm="{sm}_hifi",
        ),
    output:
        status=join(OUTPUT_DIR, "nucflag", "{sm}_status_{annot}.tsv"),
    params:
        awk_cmd=get_awk_cmd_by_annot,
        calls=lambda wc, input: (
            f"<(grep -v chrY {input.calls})" if "CHM13v2.0" in wc.sm else input.calls
        ),
    conda:
        # TODO: Update to fix bug in status
        # See https://github.com/logsdon-lab/NucFlag/pull/64
        "../Snakemake-NucFlag/workflow/env/nucflag.yaml"
    shell:
        """
        nucflag status \
            -i {params.calls} \
            -b <({params.awk_cmd} <(grep -v chrY {input.annot})) \
            -g name > {output.status}
        """


rule plot_qv_liftover_length_by_annot:
    input:
        statuses=expand(
            rules.generate_status_by_annot.output,
            annot="{annot}",
            sm=ASM_COLORS.keys(),
        ),
        group_colors=lambda wc: ANNOT_COLORS[wc.annot],
    output:
        plot=multiext(
            join(OUTPUT_DIR, "nucflag", "all_status_breakdown_by_{annot}"),
            ".png",
            ".pdf",
            ".tsv",
        ),
    params:
        script="workflow/scripts/compare_assembly/plot_status_by_annot.py",
        labels=" ".join(ASM_COLORS.keys()),
        colors=" ".join(ASM_COLORS.values()),
        output_prefix=lambda wc, output: splitext(output[0])[0],
        title=lambda wc: ANNOT_TITLE[wc.annot],
    shell:
        """
        python {params.script} \
        -i {input.statuses} \
        -l {params.labels} \
        -c {params.colors} \
        -g {input.group_colors} \
        -o {params.output_prefix} \
        --title '{params.title}'
        """


rule compare_all:
    input:
        rules.plot_compare_misassemblies.output,
        expand(
            rules.plot_qv_liftover_length_by_annot.output, annot=["segdups", "censat"]
        ),
        expand(
            rules.generate_status_by_annot.output,
            annot=["segdups", "censat"],
            sm=ASM_COLORS.keys(),
        ),
        rules.aggregate_call_num.output,
