bedtools intersect \
    -a <(grep -v good results/curated/HG002_v1.0.1_misassemblies.bed) \
    -b <(awk -v OFS="\t" '{ print $1, $2, $3, $4, $9}' /project/logsdon_shared/projects/Keith/NucFlagPaper/results/hg002v1.0.1.cenSatv2.0.noheader.bed) \
    -wb > results/curated/seq_ctx/ovl_censat.bed
bedtools intersect \
    -a <(grep -v good results/curated/HG002_v1.0.1_misassemblies.bed) \
    -b <(awk -v OFS="\t" '{ print $1, $2, $3, $37, $9}' /project/logsdon_shared/projects/Keith/NucFlagPaper/results/HG002.SDs.010624.45col.bed) \
    -wb > results/curated/seq_ctx/ovl_segdup.bed
