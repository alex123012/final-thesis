#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys
from utils import retrieve_name_from_filename, read_csv

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_ref_alt_heatmap(file_name: str, result_heatmap: str):
    name = retrieve_name_from_filename(file_name)
    df = read_csv(file_name)
    data = pd.crosstab(df.REF, df.ALT)

    fig, ax = plt.subplots()
    fig.suptitle(f"Editing events count in {name}")
    sns.heatmap(data=data, cmap="Blues", ax=ax, annot=True)
    fig.savefig(result_heatmap)


def main():
    file_name = sys.argv[1]
    result_heatmap = sys.argv[2]
    plot_ref_alt_heatmap(file_name=file_name, result_heatmap=result_heatmap)


if __name__ == "__main__":
    main()
