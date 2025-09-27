#!/bin/bash
set -euo pipefail

chrom=$1

bedtools intersect \
    -a <(grep ${chrom} /project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/HG002_v1.0.1_misassemblies.bed) \
    -b <(grep -v true_positive /project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/nucflag_v1.0.1_missed/missed_calls.tsv | grep ${chrom}) | \
    grep -v good | \
    awk -v OFS="\t" '{ if ($5 != 0 && $4 == "misjoin") { $4="indel"}; print}' | \
    cut -f 4 | \
    sort | \
    uniq -c
