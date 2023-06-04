#! /usr/bin/env python3

"""plot figures for mirna editing events"""

import sys
from typing import Dict, Callable

from utils import retrieve_name_from_filename, read_csv
from heatmap import plot_ref_alt_heatmap
from barplot import plot_ref_alt_barplot, HUE_COL, df_for_barplot
from histogram import plot_vaf_histogram

import pandas as pd


def concat_df_dict(
    df_dict: Dict[str, pd.DataFrame],
    func: Callable[[pd.DataFrame], pd.DataFrame] = None,
    hue: str = HUE_COL,
) -> pd.DataFrame:
    data = pd.DataFrame({})
    for file_name, df in df_dict.items():
        if func is not None:
            df = func(df)
        df[hue] = retrieve_name_from_filename(file_name)
        data = pd.concat([data, df])
    return data.drop_duplicates(data.columns[:-1])


def main():
    additional_figname = sys.argv[1]
    result_barplot = sys.argv[2]
    result_heatmap = sys.argv[3]
    result_hist = sys.argv[4]
    table_files = sys.argv[5:]

    sup_title = f"all samples{additional_figname}"

    df_dict = {file_name: read_csv(file_name) for file_name in table_files}

    data_heatmap = concat_df_dict(df_dict)
    plot_ref_alt_heatmap(df=data_heatmap, name=sup_title, result_heatmap=result_heatmap)

    data_histplot = concat_df_dict(df_dict=df_dict)
    plot_vaf_histogram(df=data_histplot, name=sup_title, result_hist=result_hist)

    data_barplot = concat_df_dict(df_dict=df_dict, func=df_for_barplot)
    plot_ref_alt_barplot(
        df=data_barplot,
        name=sup_title,
        result_barplot=result_barplot,
        hue=HUE_COL,
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
