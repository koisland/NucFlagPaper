#!/bin/bash

set -euo pipefail

nucflag -i /project/logsdon_shared/projects/Keith/NucFlagPaper/results/aligners/mm2/HG002_ont_r10.bam \
    -b <(printf "chr1_MATERNAL\t1\t5000000") \
    -d results/nucflag_ont_r10_w1_i1 \
    -c config/nucflag_ont_r10_w1_i1.toml

nucflag -i /project/logsdon_shared/projects/Keith/NucFlagPaper/results/aligners/mm2/HG002_ont_r10.bam \
    -b <(printf "chr1_MATERNAL\t1\t5000000") \
    -d results/nucflag_ont_r10_w11_i1 \
    -c config/nucflag_ont_r10_w11_i1.toml

nucflag -i /project/logsdon_shared/projects/Keith/NucFlagPaper/results/aligners/mm2/HG002_ont_r10.bam \
    -b <(printf "chr1_MATERNAL\t1\t5000000") \
    -d results/nucflag_ont_r10_w11_i2 \
    -c config/nucflag_ont_r10_w11_i2.toml
