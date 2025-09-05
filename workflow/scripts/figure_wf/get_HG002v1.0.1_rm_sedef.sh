#!/bin/bash

set -euo pipefail

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/segdups/HG002.SDs.010624.45col.bb -P results/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.0.1_v2.0/hg002v1.0.1.cenSatv2.0.noheader.bb -P results/
wget https://42basepairs.com/download/s3/human-pangenomics/T2T/HG002/assemblies/annotation/cytoband/cytoBand.hg002v1.0.bb -P results/
