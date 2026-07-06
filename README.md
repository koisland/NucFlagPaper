# NucFlag Paper
Code for NucFlag paper.

## Getting started
Create conda environment.
```bash
conda create --name NucFlagPaper bioconda::snakemake==9.5.0
```

To run:
```bash
snakemake -np --configfile config/config.yaml
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
* [ ] - Version and use Zenodo
* [ ] - Download CHM13 data automatically for assembly comparison.
* [ ] - Refactor to run in parts.
* [ ] - Rename so better structured (ex. `1-*`, `2-*` )
* [ ] - Test using pixi
