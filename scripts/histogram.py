#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import retrieve_name_from_filename, read_csv

VAF_COL = "VAF"


def plot_vaf_histogram(
    df: pd.DataFrame, name: str, result_hist: str, x: str = VAF_COL, hue: str = None
):
    data = df.copy()
    if data.shape[0] < 4:
        return

    sns.set_context("poster")
    fig, ax = plt.subplots(figsize=(20, 14))
    fig.suptitle(f"VAF distribution in {name}")
    sns.histplot(
        data=data,
        x=x,
        ax=ax,
        hue=hue if hue else None,
        kde=True,
    )
    fig.savefig(result_hist)


def main():
    file_name = sys.argv[1]
    result_hist = sys.argv[2]

    name = retrieve_name_from_filename(file_name)
    df = read_csv(file_name)
    plot_vaf_histogram(df=df, name=name, x=VAF_COL, result_hist=result_hist)


if __name__ == "__main__":
    main()
