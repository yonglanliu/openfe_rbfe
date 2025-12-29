"""
Docstring for openfe_app.json_analysis
created by: Yonglan Liu
Date: 2025/12/22
"""

import json
import numpy as np
import pandas as pd

# --------------------------------------------------
# Function to extract deltaG: all repeats and final
# --------------------------------------------------
def extract_dg(data):
    repeats = list(data["unit_results"].keys())
    unit = []
    dg = []
    uncertainty =[]
    source_keys = []
    for i, rep in enumerate(repeats):
        source_key = data["unit_results"][rep]["source_key"]
        rep_result = data["unit_results"][rep]["outputs"]
        rep_dg = rep_result["unit_estimate"]["magnitude"]
        rep_uncertainty = rep_result["unit_estimate_error"]["magnitude"]
        unit.append(f'rep{i}')
        dg.append(rep_dg)
        uncertainty.append(rep_uncertainty)
        source_keys.append(source_key)
    #print(np.average(dg))
    #print(np.average(uncertainty))

    unit.append("average_of_replicates")
    dg.append(data["estimate"]["magnitude"])
    uncertainty.append(data["uncertainty"]["magnitude"])
    source_keys.append("Final")

    df_final = pd.DataFrame({
    "Unit": unit,
    "DG (kcal/mol)": dg,
    "uncertainty (kcal/mol)": uncertainty,
    "source_key":source_keys
    })
    return df_final

# ----------------------------------------
#  Function to decode JSON ndarray data
# ----------------------------------------
def decode_openfe_ndarray(obj):
    """
    Decode OpenFE/gufe JSON ndarray objects.
    Handles zstd-compressed payloads (magic 28 B5 2F FD) and raw buffers.
    """
    if not (isinstance(obj, dict) and obj.get("__class__") == "ndarray"):
        return obj

    dtype = np.dtype(obj["dtype"])
    shape = tuple(obj["shape"])
    expected_n = int(np.prod(shape))
    expected_bytes = expected_n * dtype.itemsize

    b = obj["bytes"]["latin-1"].encode("latin-1")

    # Zstandard magic number
    is_zstd = len(b) >= 4 and b[:4] == b"\x28\xb5\x2f\xfd"
    if is_zstd:
        import zstandard as zstd
        dctx = zstd.ZstdDecompressor()
        raw = dctx.decompress(b)

        if len(raw) != expected_bytes:
            raise ValueError(f"Decompressed size {len(raw)} != expected {expected_bytes}")

        return np.frombuffer(raw, dtype=dtype).reshape(shape)

    # If not zstd, assume raw buffer (must match expected size)
    if len(b) != expected_bytes:
        raise ValueError(f"Raw buffer size {len(b)} != expected {expected_bytes}")

    return np.frombuffer(b, dtype=dtype).reshape(shape)

# ---------------------------------------------------------------------
# Function to extract forward/reverse convergence data and eigenvalues
# ---------------------------------------------------------------------
def extract_f_and_r_convergence(rep):
    #print((rep["outputs"]))
    fre = rep["outputs"]["forward_and_reverse_energies"]
    fractions = decode_openfe_ndarray(fre["fractions"])  # usually already a plain list

    forward_DG  = decode_openfe_ndarray(fre["forward_DGs"]["magnitude"])
    reverse_DG  = decode_openfe_ndarray(fre["reverse_DGs"]["magnitude"])
    forward_err = decode_openfe_ndarray(fre["forward_dDGs"]["magnitude"])
    reverse_err = decode_openfe_ndarray(fre["reverse_dDGs"]["magnitude"])

    df_dg = pd.DataFrame({
        "fraction": fractions,
        "forward_DG (kcal/mol)": forward_DG,
        "reverse_DG (kcal/mol)": reverse_DG,
        "forward_error (kcal/mol)": forward_err,
        "reverse_error (kcal/mol)": reverse_err,
    })
    return df_dg

# -----------------------------------------------------------------
# Function to extract replica exchange matrix data and eigenvalues
# -----------------------------------------------------------------
def extract_rex(rep):
    rex = rep["outputs"]["replica_exchange_statistics"]
    
    rex_matrix = decode_openfe_ndarray(rex["matrix"])
    rex_eig = decode_openfe_ndarray(rex["eigenvalues"])

    n_lambda = rex_matrix.shape[0]
    """
    df_rex = pd.DataFrame(
        rex_matrix,
        index=[f"λ{i}" for i in range(n_lambda)],
        columns=[f"λ{i}" for i in range(n_lambda)],
    )
    """
    df_rex_eig = pd.DataFrame({
    "eigenvalue": rex_eig,
    "timescale_index": np.arange(len(rex_eig))
    })
    return np.array(rex_matrix, dtype=float), df_rex_eig

# ------------------------------------------------------
#  Function to extract mbar matrix and and eigenvalues
# ------------------------------------------------------
def extract_mbar(rep):
    mbar = rep["outputs"]["unit_mbar_overlap"]
    mbar_matrix = decode_openfe_ndarray(mbar["matrix"])
    mbar_eig = decode_openfe_ndarray(mbar["eigenvalues"])
    n_lambda = mbar_matrix.shape[0]
    """
    df_mbar = pd.DataFrame(
        mbar_matrix,
        index=[f"λ{i}" for i in range(n_lambda)],
        columns=[f"λ{i}" for i in range(n_lambda)],
    )
    """
    df_mbar_eig = pd.DataFrame({
    "eigenvalue": mbar_eig,
    "index": np.arange(len(mbar_eig))
    })
    return np.array(mbar_matrix), df_mbar_eig

"""
def main():
    json_file = "openfe_run_by_me/solvent/results_solvent.json"
    # Load json file
    with open(json_file) as f:
        data = json.load(f)
    # Json file contains keys: ['estimate', 'uncertainty', 'protocol_result', 'unit_results']

    # Extract data
    protocol_result = data["protocol_result"]
    repeats = protocol_result["data"]
    repeat_ids = list(repeats.keys())
    print(repeat_ids)

    dg_final = extract_dg(data)

"""