import os
import pandas as pd


def retrieve_name_from_filename(file_name: str) -> str:
    bname = os.path.basename(file_name)
    name = os.path.splitext(bname)[0]
    return name.lstrip("result-")


def read_csv(file_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_name, sep="\t")
    df.ALT = df.ALT.str.split(",")
    df = df.explode("ALT")
    return df
