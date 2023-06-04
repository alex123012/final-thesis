#! /usr/bin/env python3

"""plot figures for mirna editing events"""


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from barplot import X_COL, Y_COL


EDITED_IN_COL = "EDITED IN"


def count_seed_editing_events(df: pd.DataFrame) -> pd.DataFrame:
    if df.shape[0] < 1:
        return pd.DataFrame({})

    seed_edit = set()
    for _, mirline in df[df.TYPE == "miRNA"].iterrows():
        seq_len = len(mirline.REF_SEQ)
        edit_pos = mirline.REF_SEQ_POS + 1

        pass_strand_seed_pos = (
            mirline.STRAND == "-" and edit_pos < seq_len and edit_pos > seq_len - 8
        )
        prim_strand_seed_pos = mirline.STRAND == "+" and edit_pos > 1 and edit_pos < 8

        if pass_strand_seed_pos or prim_strand_seed_pos:
            seed_edit.add(mirline.POS)

    df[EDITED_IN_COL] = "outside seed region"
    for mir_pos in seed_edit:
        df[EDITED_IN_COL] = df.apply(
            in_seed_region(mir_pos),
            axis=1,
        )

    df[X_COL] = df.apply(lambda row: row.REF + " -> " + row.ALT, axis=1)
    return df


def in_seed_region(pos: int) -> str:
    def func(row: pd.Series) -> str:
        return (
            "seed region" if pos > row.START and pos < row.END else row[EDITED_IN_COL]
        )

    return func


def plot_seed_editing_events_count(df: pd.DataFrame, name: str, result_image: str):
    data = (
        pd.crosstab(df[EDITED_IN_COL], df[X_COL])
        .melt(ignore_index=False, value_name=Y_COL)
        .reset_index()
    )

    sns.set_context("poster")
    fig, ax = plt.subplots(figsize=(20, 14))
    fig.suptitle(f"Editing events in seed regions and outside seed regions in\n{name}")

    palette = dict(
        zip(
            ["seed region", "outside seed region"],
            sns.color_palette(n_colors=2),
        )
    )
    sns.barplot(
        data=data.sort_values(Y_COL),
        x=X_COL,
        y=Y_COL,
        hue=EDITED_IN_COL,
        ax=ax,
        palette=palette,
    )
    fig.savefig(result_image)
