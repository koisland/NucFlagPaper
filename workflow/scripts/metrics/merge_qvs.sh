#!/bin/bash

set -euo pipefail

awk -v OFS="\t" 'FNR > 1 {
    match(FILENAME, "(v[^_/]+)_nucflag", vs);
    print $1, $2, $3, $4, vs[1], "NucFlag v1.0", "PacBio HiFi"
}' /project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/qv/*_qv.bed
awk -v OFS="\t" '{
    match(FILENAME, "qv/HG002_(.+?)/", arr);
    print $1, $2, $3, $4, arr[1], "Merqury v1.4", "Element"
}' /project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/qv/HG002_*/HG002_v*.v*.qv
