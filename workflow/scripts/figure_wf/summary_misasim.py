import polars as pl


df_nucflag = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/misasim/summary/nucflag_HG002_*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)
df_flagger = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/misasim/summary/flagger_HG002_*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)
df_inspector = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/misasim/summary/inspector_HG002_*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)


def calculate_metrics(df: pl.DataFrame) -> tuple[float, float, float]:
    prec = df["precision"].mean()
    recall = df["recall"].mean()
    f1 = (2 * (prec * recall)) / (prec + recall)
    return prec, recall, f1


nf_prec, nf_recall, nf_f1 = calculate_metrics(df_nucflag)
flg_prec, flg_recall, flg_f1 = calculate_metrics(df_flagger)
insp_prec, insp_recall, insp_f1 = calculate_metrics(df_inspector)

print(nf_prec, nf_recall, nf_f1)
print(flg_prec, flg_recall, flg_f1)
print(insp_prec, insp_recall, insp_f1)
