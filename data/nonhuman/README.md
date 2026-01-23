# Centromere locations
Get centromere coordinates in Col-XJTU paper
* `/project/logsdon_shared/projects/Keith/NucFlagPaper/data/nonhuman/Col-XJTU_centromeres.bed`
* In section, *"Assembly validation of CEN180 arrays"*

Get centromere coordinates in Col-PEK paper.
* `/project/logsdon_shared/projects/Keith/NucFlagPaper/data/nonhuman/Col-PEK_centromeres.bed`
* In section, *"Supplemental Table 4. Merqury evaluation of the regions spanned by the five CEN180 arrays and unmapped ONT raw reads of Col-PEK"*

Get centromere coordinates in Col-CEN paper.
* `/project/logsdon_shared/projects/Keith/NucFlagPaper/data/nonhuman/Col-CEN_centromeres.bed`
* In section, *"Table S1. Consensus quality (QV) score of the Col-CEN Arabidopsis genome assembly"*

## Liftover
Alternatively, use `asm-to-ref` pipeline to align Col-PEK to Col-XJTU. Done in paper workflow.

Invert PAF of Col-XJTU (query) to Col-PEK (ref) with `rustybam`.
```bash
rb invert /project/logsdon_shared/projects/Keith/NucFlagPaper/results/nonhuman/saffire/Col-PEK/paf/Col-C.paf > /project/logsdon_shared/projects/Keith/NucFlagPaper/exp/nonhuman/Col-XJTU.paf
```

Liftover coordinates to Col-PEK with `impg`.
```bash
impg query \
-p /project/logsdon_shared/projects/Keith/NucFlagPaper/exp/nonhuman/Col-XJTU.paf \
-b /project/logsdon_shared/projects/Keith/NucFlagPaper/data/nonhuman/Col-XJTU_centromeres.bed | \
bedtools groupby -i - -g 1 -c 2,3 -o min,max > data/nonhuman/Col-PEK_centromeres.bed
```

Then liftover Col-PEK to Col-CEN with `impg`.
```bash
impg query \
-p /project/logsdon_shared/projects/Keith/NucFlagPaper/results/nonhuman/saffire/Col-PEK/paf/Col-CEN.paf \
-b data/nonhuman/Col-PEK_centromeres.bed | \
bedtools groupby -i - -g 1 -c 2,3 -o min,max > data/nonhuman/Col-CEN_centromeres.bed
```
