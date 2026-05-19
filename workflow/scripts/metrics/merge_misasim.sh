#!/bin/bash

set -euo pipefail

indir=/project/logsdon_shared/projects/Keith/NucFlagPaper/results/misasim/summary

awk -v OFS="\t" 'FNR > 1 {
    match(FILENAME, "([^/]*?).tsv", arr);
    split(arr[1], fin_arr, "_");
    print fin_arr[1], $0
}' ${indir}/*_hifi.tsv ${indir}/*_ont_r10.tsv | \
sed -e 's|flagger|HMM-Flagger v1.1.0|g' \
    -e 's|inspector|Inspector v1.3|g' \
    -e 's|nucflag|NucFlag v1.0.0|g' \
    -e 's|hifi|PacBio HiFi|g' \
    -e 's|ont_r10|ONT R10|g'
