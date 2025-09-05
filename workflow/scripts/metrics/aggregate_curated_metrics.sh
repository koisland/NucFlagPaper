#!/bin/bash

set -euo pipefail

awk -v OFS="\t" \
    'FNR > 1 {
        match(FILENAME, "/(nucflag|flagger|inspector)_v(.+).tsv", mtch);
        print mtch[1], mtch[2], $0
    }' /project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/*.tsv
