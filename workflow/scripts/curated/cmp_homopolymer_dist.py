import argparse
import polars as pl
from scipy.stats import ks_2samp, false_discovery_control


def df_counts(df: pl.DataFrame) -> pl.DataFrame:
    return df.group_by(["length"]).agg(pl.col("count").sum()).sort("length")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--grp_a", nargs="+")
    ap.add_argument("-b", "--grp_b", nargs="+")
    ap.add_argument("-l", "--labels", nargs="+")
    args = ap.parse_args()

    rows = []
    for lbl, a, b in zip(args.labels, args.grp_a, args.grp_b):
        df_a = pl.read_csv(
            a, separator="\t", has_header=False, new_columns=["nt", "length", "count"]
        ).filter(pl.col("length") > 0)
        df_b = pl.read_csv(
            b, separator="\t", has_header=False, new_columns=["nt", "length", "count"]
        ).filter(pl.col("length") > 0)

        for nt in ("A", "T", "G", "C"):
            df_nt_a = df_a.filter(pl.col("nt").eq(pl.lit(nt)))
            df_nt_b = df_b.filter(pl.col("nt").eq(pl.lit(nt)))
            df_counts_a = df_counts(df_nt_a)
            df_counts_b = df_counts(df_nt_b)
            # Use kolmogorov-smirnov to compare if distributions are identical.
            res = ks_2samp(df_counts_a["count"], df_counts_b["count"])
            rows.append((lbl, nt, res.pvalue, res.statistic))

        df_counts_a_all = df_counts(df_a)
        df_counts_b_all = df_counts(df_b)
        res_all = ks_2samp(df_counts_a_all["count"], df_counts_b_all["count"])

        rows.append((lbl, "all", res_all.pvalue, res_all.statistic))

    lbls, nts, pvals, stats = zip(*rows)
    # FDR control
    pvals_new = false_discovery_control(pvals)
    for lbl, nt, pval, stat in zip(lbls, nts, pvals_new, stats):
        print(lbl, nt, pval, stat, sep="\t")


if __name__ == "__main__":
    raise SystemExit(main())
