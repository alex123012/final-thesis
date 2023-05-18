#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys
from typing import Dict, Callable

from utils import retrieve_name_from_filename, read_csv
from heatmap import plot_ref_alt_heatmap
from barplot import plot_ref_alt_barplot, HUE_COL, df_for_barplot

import pandas as pd


def concat_df_dict(
    df_dict: Dict[str, pd.DataFrame],
    func: Callable[[pd.DataFrame], pd.DataFrame] = None,
) -> pd.DataFrame:
    data = pd.DataFrame({})
    for file_name, df in df_dict.items():
        if func is not None:
            df = func(df)
        df[HUE_COL] = retrieve_name_from_filename(file_name)
        data = pd.concat([data, df])
    return data.drop_duplicates(data.columns[:-1])


def main():
    result_barplot = sys.argv[1]
    result_heatmap = sys.argv[2]
    table_files = sys.argv[3:]

    df_dict = {file_name: read_csv(file_name) for file_name in table_files}

    data_heatmap = concat_df_dict(df_dict)
    plot_ref_alt_heatmap(
        df=data_heatmap, name="All samples", result_heatmap=result_heatmap
    )

    data_barplot = concat_df_dict(df_dict=df_dict, func=df_for_barplot)
    plot_ref_alt_barplot(
        df=data_barplot, name="All samples", result_barplot=result_barplot, hue=HUE_COL
    )

    print("Events found: ", data_heatmap.shape[0])
    print("VAF mean: ", data_heatmap.VAF.mean())
    print("VAF median: ", data_heatmap.VAF.median())
    print(
        "VAF mean without canonical: ",
        data_heatmap[
            ((data_heatmap.ALT != "U") & (data_heatmap.REF != "C"))
            & ((data_heatmap.ALT != "G") & (data_heatmap.REF != "A"))
        ].VAF.mean(),
    )
    print(
        "VAF median without canonical: ",
        data_heatmap[
            ((data_heatmap.ALT != "U") & (data_heatmap.REF != "C"))
            & ((data_heatmap.ALT != "G") & (data_heatmap.REF != "A"))
        ].VAF.median(),
    )


if __name__ == "__main__":
    main()
