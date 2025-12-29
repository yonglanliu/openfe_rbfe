from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def plot_dg_table(df: pd.DataFrame):
    """
    For `openfe gather --report dg`, the table typically includes a ligand identifier,
    DG estimate, and uncertainty columns (names can vary slightly by version).
    We'll try to infer likely columns.
    """
    # Try common column name guesses
    cols = [c.lower() for c in df.columns]
    df2 = df.copy()
    df2.columns = cols

    # pick x label column
    label_col = None
    for cand in ["ligand", "name", "ligand_name", "mol", "molecule"]:
        if cand in df2.columns:
            label_col = cand
            break
    if label_col is None:
        label_col = df2.columns[0]

    # pick dg column
    dg_col = None
    for cand in ["dg", "delta_g", "g", "free_energy"]:
        if cand in df2.columns:
            dg_col = cand
            break
    # fallback: first numeric col
    if dg_col is None:
        for c in df2.columns:
            if pd.api.types.is_numeric_dtype(df2[c]):
                dg_col = c
                break

    # pick uncertainty column
    err_col = None
    for cand in ["uncertainty", "dg_std", "sigma", "error", "stderr"]:
        if cand in df2.columns:
            err_col = cand
            break

    # Plot
    fig, ax = plt.subplots()
    x = np.arange(len(df2))
    y = df2[dg_col].astype(float).values
    if err_col and pd.api.types.is_numeric_dtype(df2[err_col]):
        yerr = df2[err_col].astype(float).values
        ax.errorbar(x, y, yerr=yerr, fmt="o")
    else:
        ax.plot(x, y, "o")

    ax.set_xticks(x)
    ax.set_xticklabels(df2[label_col].astype(str).values, rotation=90)
    ax.set_ylabel(dg_col)
    ax.set_title("Gathered DG results")
    fig.tight_layout()
    return fig
