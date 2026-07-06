# NucFlag Paper
Code for NucFlag paper and benchmarking.

## Getting started
Clone repo and its submodules.
```bash
git clone https://github.com/koisland/NucFlagPaper.git --recursive
```

Create conda environment.
```bash
conda create --name NucFlagPaper bioconda::snakemake==9.5.0
```

## Usage
To run all analyses for paper:
```bash
which apptainer
# Replace profiles/default/ with preferred execution method and rop -e.
snakemake -np --configfile config/config.yaml -e local
```

> [!NOTE]
> This requires a lot of alignment jobs and disk space. Prefer the options below.

<details>
<summary>Jobs run</summary>

```
Job stats:
job                                                       count
------------------------------------------------------  -------
all                                                           1
aln_align_reads_to_asm                                        6
aln_cmp_calculate_precision_recall                            6
aln_cmp_nf_mm2_aln_align_reads_to_asm                       374
aln_cmp_nf_mm2_aln_merge_asm_files                           44
aln_cmp_nf_mm2_aln_merge_read_asm_alignments                 44
aln_cmp_nf_mm2_check_asm_nucflag                             44
aln_cmp_nf_pbmm2_aln_align_reads_to_asm                     374
aln_cmp_nf_pbmm2_aln_merge_asm_files                         44
aln_cmp_nf_pbmm2_aln_merge_read_asm_alignments               44
aln_cmp_nf_pbmm2_check_asm_nucflag                           44
aln_cmp_nf_wm2_aln_align_reads_to_asm                       374
aln_cmp_nf_wm2_aln_get_repetitive_kmers                      44
aln_cmp_nf_wm2_aln_merge_asm_files                           44
aln_cmp_nf_wm2_aln_merge_read_asm_alignments                 44
aln_cmp_nf_wm2_check_asm_nucflag                             44
aln_cmp_plot_precision_recall                                 1
aln_merge_asm_files                                           6
aln_merge_read_asm_alignments                                 6
bam_to_cov                                                  132
check_asm_nucflag                                             6
cmp_asm_aggregate_call_num                                    1
cmp_asm_asm_ref_SafFire                                       3
cmp_asm_asm_ref_alignment                                     3
cmp_asm_asm_ref_bam_to_paf                                    3
cmp_asm_asm_ref_trim_and_break_paf                            3
cmp_asm_calculate_stats                                       1
cmp_asm_create_hor_array_bed                                  3
cmp_asm_denovo_cleanup_tmp_fastq                              3
cmp_asm_denovo_generate_dtype_fofn                            6
cmp_asm_denovo_generate_original_dtype_fofn                   6
cmp_asm_denovo_generate_summary_stats                         3
cmp_asm_denovo_hifiasm_output                                 1
cmp_asm_denovo_run_assembler                                  3
cmp_asm_denovo_verkko_output                                  2
cmp_asm_download_annotations                                  1
cmp_asm_draw_ng50                                             1
cmp_asm_find_all_lcr_regions                                  4
cmp_asm_generate_filtered_calls                               4
cmp_asm_generate_status_bed                                   3
cmp_asm_generate_status_by_annot                              8
cmp_asm_intersect_bubbles_w_calls                             4
cmp_asm_label_lcr_regions_and_pos                             4
cmp_asm_liftover_annotations                                  3
cmp_asm_merge_beds                                            1
cmp_asm_nf_denovo_aln_align_reads_to_asm                     20
cmp_asm_nf_denovo_aln_merge_asm_files                         8
cmp_asm_nf_denovo_aln_merge_read_asm_alignments               8
cmp_asm_nf_denovo_check_asm_nucflag                           8
cmp_asm_nf_denovo_generate_breakdown_plot                     8
cmp_asm_nf_denovo_generate_ideogram                           8
cmp_asm_overlap_cmp_dist                                      1
cmp_asm_paf2chain                                             3
cmp_asm_plot_compare_misassemblies                            1
cmp_asm_plot_heatmap_chm13_impg                               1
cmp_asm_plot_qv_liftover_length_by_annot                      2
cmp_asm_query_chm13_impg                                      1
cmp_asm_rename_assemblies                                     4
cmp_asm_run_cenmap                                            3
cmp_asm_vg_ava_asm                                            1
cmp_asm_vg_generate_variation_graph_output                    1
cmp_asm_vg_merge_asm                                          1
cmp_asm_vg_split_regions                                      1
create_wg_bed                                               132
download_all                                                  1
download_files                                               25
hg002_curated_aggregate_call_num                              1
hg002_curated_align_element                                   3
hg002_curated_alignment_element_remove_dup                    3
hg002_curated_calculate_nucflag_qv                            3
hg002_curated_calculate_precision_recall                     12
hg002_curated_compare_nf_genome_dist                          1
hg002_curated_convert_annotations_to_bed                      1
hg002_curated_convert_vcf_to_bed                              3
hg002_curated_download_annotation_data                        1
hg002_curated_download_curated_asm                            3
hg002_curated_download_element                                1
hg002_curated_download_vcf                                    3
hg002_curated_find_all_lcr_regions                            3
hg002_curated_generate_index                                  3
hg002_curated_get_consensus_nucflag                           3
hg002_curated_label_lcr_regions_and_pos                       3
hg002_curated_liftover_cytobands_from_chm13                   1
hg002_curated_liftover_from_chm13_alignment                   1
hg002_curated_liftover_from_chm13_bam_to_paf                  1
hg002_curated_liftover_from_chm13_paf2chain                   1
hg002_curated_merge_element_mapq_bigwigs                      3
hg002_curated_overlap_annotation_data                         4
hg002_curated_plot_call_venn                                  3
hg002_curated_plot_ideogram                                   2
hg002_curated_plot_ideogram_chrom                             6
hg002_curated_plot_mapq_with_consensus                        3
hg002_curated_plot_nucflag_merqury_qv_plot                    1
hg002_curated_plot_qv_functions                               1
hg002_curated_plot_seq_context                                4
hg002_curated_plot_stats_f1_all                               1
hg002_curated_qv_count_kmers                                  6
hg002_curated_qv_estimate_qv                                  3
hg002_curated_qv_merge_dbs                                    3
hg002_curated_recalculate_precision_recall_w_intersect        3
hg002_curated_versioned_aln_align_reads_to_asm                9
hg002_curated_versioned_aln_merge_asm_files                   3
hg002_curated_versioned_aln_merge_read_asm_alignments         3
hg002_curated_versioned_bam_to_cov                            3
hg002_curated_versioned_check_asm_nucflag                     3
hg002_curated_versioned_create_wg_bed                         3
hg002_curated_versioned_element_check_asm_nucflag             3
hg002_curated_versioned_filter_vcf_calls                      3
hg002_curated_versioned_generate_breakdown_plot               3
hg002_curated_versioned_generate_ideogram                     3
hg002_curated_versioned_make_annot_json                       3
hg002_curated_versioned_merge_calls                           3
hg002_curated_versioned_run_deepvariant                       3
hg002_curated_versioned_run_flagger                           3
hg002_curated_versioned_run_inspector                         3
hprc_Rhodonite_DupMasker                                     80
hprc_Rhodonite_run_DupMasker_step_2                         160
hprc_Rhodonite_run_DupMasker_step_3                         160
hprc_Rhodonite_run_split_RepeatMasker                       160
hprc_Rhodonite_setup_DupMasker                                1
hprc_Rhodonite_setup_RepeatMasker                             1
hprc_Rhodonite_split_fasta                                   80
hprc_asm_ref_hprc_r1_r2_SafFire                              80
hprc_asm_ref_hprc_r1_r2_alignment                            80
hprc_asm_ref_hprc_r1_r2_aln_to_bed                           80
hprc_asm_ref_hprc_r1_r2_bam_to_paf                           80
hprc_asm_ref_hprc_r1_r2_trim_and_break_paf                   80
hprc_calculate_qv                                            80
hprc_calculate_release_ng50                                   2
hprc_download_aligned_bams                                   80
hprc_download_assemblies                                     80
hprc_download_r2_cens                                        40
hprc_generate_cens_bed                                       40
hprc_intersect_dupmasker_nucflag                             80
hprc_merge_intervals_preserve_ort                            80
hprc_merge_statuses_and_label                                 2
hprc_plot_err_comparison                                      1
hprc_plot_qv_comparison                                       2
hprc_plot_r1_r2_chm13_coords                                  1
hprc_query_w_impg_sm_to_chm13_paf                            80
hprc_run_nucflag                                             80
hprc_status_count_by_region                                  80
hprc_subset_smn                                              80
make_annot_json                                             132
merge_calls                                                 132
nonhuman_agg_calls_cen                                        2
nonhuman_agg_qv_stats_against_reported                        1
nonhuman_asm_ref_SafFire                                      2
nonhuman_asm_ref_alignment                                    2
nonhuman_asm_ref_bam_to_paf                                   2
nonhuman_asm_ref_trim_and_break_paf                           2
nonhuman_calculate_qv                                         3
nonhuman_create_groups                                        3
nonhuman_download_data                                        9
nonhuman_generate_ideogram_all                                1
nonhuman_plot_agg_calls_cen                                   1
nonhuman_status_groups                                        6
run_flagger                                                 132
run_inspector                                               132
sim_calculate_precision_recall                                2
sim_calculate_precision_recall_flagger                        2
sim_calculate_precision_recall_inspector                      2
sim_create_tmp_bam                                          132
sim_downsample_reads                                          4
sim_generate_benchmarks                                       1
sim_ms_aln_align_reads_to_asm                               462
sim_ms_aln_merge_asm_files                                  132
sim_ms_aln_merge_read_asm_alignments                        132
sim_ms_check_asm_nucflag                                    132
sim_ms_compile_misasim                                        1
sim_ms_generate_inv_bed                                      42
sim_ms_generate_misassemblies                                42
sim_ms_index_misassemblies_fa                                42
sim_ms_split_asm_misasim                                      2
sim_plot_precision_recall                                     2
sim_plot_precision_recall_xy                                  2
sim_plot_stats_f1_all                                         2
total                                                      5731
```
</details>

To run just assembly error simulation.
```bash
TODO
```

To run assembly benchmarking for hifiasm and verkko with different versions.
```bash
TODO
```

## Organization
This repository contains all code used in the NucFlag paper. Its organized as a standard Snakemake workflow.
```
workflow/
├── envs
│   ├── cenmap.yaml
│   ├── compare_assembly.yaml
│   ├── curated.yaml
│   ├── download.yaml
│   ├── inspector.yaml
│   └── tools.yaml
├── profiles
│   ├── default
│   └── lpc_all
├── rules
│   ├── aligners
│   ├── asm-to-reference-alignment
│   ├── check_curated.smk
│   ├── common.smk
│   ├── compare_aligners.smk
│   ├── compare_assembly
│   ├── compare_assembly.smk
│   ├── compare_hprc
│   ├── compare_hprc.smk
│   ├── curated
│   ├── misasim
│   ├── nonhuman
│   ├── nonhuman.smk
│   ├── simulate_misassemblies.smk
│   ├── Snakemake-NucFlag
│   └── utils
├── scripts
│   ├── compare_assembly
│   ├── compare_hprc
│   ├── curated
│   ├── figure_wf
│   ├── flagger
│   ├── Inspector
│   ├── metrics
│   ├── misasim
│   └── nonhuman
└── Snakefile
```

### Rules
Rules are split as follows. Each rule is an entrypoint with an associated subdirectory (`workflow/rules/{rule}`) and scripts directory (`workflow/scripts/{rule}`).
* `simulate_misassemblies.smk`
    * Workflow and code used to generate Figure 2
    * The subdir `workflow/rules/misasim` corresponds to this rule.
* `check_curated.smk`
    * Workflow and code used to generate Figure 3
* `compare_assembly.smk`
    * Workflow and code used to generate Figure 4
    * This workflow requires downloading CHM13 hifi data separately. This is a TODO.
* `nonhuman.smk`
    * Workflow and code used to generate Figure 5
* `compare_hprc.smk`
    * Workflow and code used to generate Figure 6
* `common.smk`
    * Helper code used for `simulate_misassemblies.smk`
* `compare_aligners.smk`
    * Workflow and code used to generate Extended Data Figure X
    * The subdir `workflow/rules/aligners` corresponds to this rule.

## Cite

## TODO
* [ ] Version and use Zenodo
* [x] Download CHM13 data automatically for assembly comparison.
* [ ] Refactor to run in parts.
* [ ] Rename so better structured (ex. `1-*`, `2-*` )
* [ ] Test using pixi
* [ ] Fork repo into separate workflow for benchmarking.
