#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys
from utils import retrieve_name_from_filename, read_csv

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_ref_alt_heatmap(df: pd.DataFrame, name: str, result_heatmap: str):
    data = pd.crosstab(df.REF, df.ALT).astype(int)
    if data.shape[0] == 0:
        return

    fig, ax = plt.subplots()
    fig.suptitle(f"Editing events count in {name}")
    sns.heatmap(data=data, cmap="Blues", ax=ax, annot=True, fmt='g')
    fig.savefig(result_heatmap)


def main():
    file_name = sys.argv[1]
    result_heatmap = sys.argv[2]

    name = retrieve_name_from_filename(file_name)
    df = read_csv(file_name)
    plot_ref_alt_heatmap(df=df, name=name, result_heatmap=result_heatmap)


if __name__ == "__main__":
    main()
