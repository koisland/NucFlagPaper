import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns


df_nucflag = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/nucflag_v*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)
df_flagger = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/flagger_v*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)
df_inspector = pl.read_csv(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/inspector_v*.tsv",
    separator="\t",
    glob=True,
    has_header=True,
)
# https://www.statology.org/what-is-a-good-f1-score/
# F1 Score = 2 * (Precision * Recall) / (Precision + Recall)


def calculate_metrics(df: pl.DataFrame) -> tuple[float, float, float]:
    prec = df["precision"].mean()
    recall = df["recall"].mean()
    f1 = 2 * (prec * recall) / (prec + recall)
    return prec, recall, f1


nf_prec, nf_recall, nf_f1 = calculate_metrics(df_nucflag)
flg_prec, flg_recall, flg_f1 = calculate_metrics(df_flagger)
insp_prec, insp_recall, insp_f1 = calculate_metrics(df_inspector)

df = pl.DataFrame(
    {
        "tool": ["HMM-Flagger v1.1.0", "Inspector v1.3", "NucFlag v1.0"] * 2,
        "percent": [flg_recall, insp_recall, nf_recall, flg_prec, insp_prec, nf_prec],
        "type": ["recall"] * 3 + ["precision"] * 3,
    }
)

ax = sns.barplot(
    df,
    x="tool",
    y="percent",
    hue="type",
    palette="pastel",
)
for bar in ax.containers:
    ax.bar_label(bar, fontsize=10, fmt=lambda x: f"{x * 100:.1f}%")

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.set_xlabel(None)

plt.savefig(
    "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/summary_curated.png"
)

print(nf_prec, nf_recall, nf_f1)
print(flg_prec, flg_recall, flg_f1)
print(insp_prec, insp_recall, insp_f1)
