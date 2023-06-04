#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys
from utils import retrieve_name_from_filename, read_csv

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


Y_COL = "Count"
X_COL = "REF -> ALT"
HUE_COL = "Samples"


def df_for_barplot(df: pd.DataFrame) -> pd.DataFrame:
    return (
        pd.crosstab(df.REF, df.ALT)
        .melt(ignore_index=False, value_name=Y_COL)
        .reset_index()
    )


def plot_ref_alt_barplot(
    df: pd.DataFrame,
    result_barplot: str,
    name: str,
    x: str = X_COL,
    y: str = Y_COL,
    hue: str = None,
):
    data = df[~(df.ALT == df.REF)].astype({Y_COL: int})
    if data.shape[0] == 0:
        return

    data[x] = data.apply(lambda row: row.REF + " -> " + row.ALT, axis=1)

    sns.set_context("poster")
    fig, ax = plt.subplots(figsize=(30, 20))
    fig.suptitle(f"Editing events count in {name}")
    order = (
        data.groupby([X_COL, Y_COL])
        .apply(lambda df: df[Y_COL].sum())
        .reset_index()
        .groupby(X_COL)
        .apply(lambda df: df[0].sum())
        .sort_values()
        .index
    )
    sns.barplot(
        data=data.sort_values(Y_COL),
        x=x,
        y=y,
        hue=hue if hue else None,
        ax=ax,
        order=order,
    )
    mx, mn = data[Y_COL].max(), data[Y_COL].min() - 5
    ax.set_yticks(
        list(
            range(
                (mn > 0 and mn - 5) or 0,
                (mx > 5 and mx + 5) or mx + 1,
                (mx > 5 and 5) or 1,
            )
        )
    )
    ax.legend(loc="upper left")
    fig.savefig(result_barplot)


def main():
    file_name = sys.argv[1]
    result_barplot = sys.argv[2]

    name = retrieve_name_from_filename(file_name)
    df = df_for_barplot(read_csv(file_name))
    plot_ref_alt_barplot(
        df=df, name=name, result_barplot=result_barplot, x=X_COL, y=Y_COL
    )


if __name__ == "__main__":
    main()
