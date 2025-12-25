import argparse

import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-x", "--qv_x", type=str, help="QV for x-axis.")
    ap.add_argument("-y", "--qv_y", type=str, help="QV for y-axis.")
    ap.add_argument("-o", "--output_prefix")

    args = ap.parse_args()
    df_qv_x = pl.read_csv(
        args.qv_x, separator="\t", has_header=False, comment_prefix="#"
    )
    df_qv_y = pl.read_csv(
        args.qv_y, separator="\t", has_header=False, comment_prefix="#"
    )
    print(df_qv_x, df_qv_y)


if __name__ == "__main__":
    raise SystemExit(main())
