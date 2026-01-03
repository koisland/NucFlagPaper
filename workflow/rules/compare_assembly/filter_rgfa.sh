#!/bin/bash

set -euo pipefail

mkdir -p results/assembly/vg/graphs_nonempty && \
find results/assembly/vg/graphs -name "*.gfa" -not -empty | \
xargs -I {} bash -c 'if [ $(grep -Po "(?<=SN:Z:)(.*?)(?=\t)" {} | sort -u | wc -l) -gt 1 ]; then cp {} results/assembly/vg/graphs_nonempty; fi;'
