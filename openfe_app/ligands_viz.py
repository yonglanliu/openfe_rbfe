from __future__ import annotations
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd
from rdkit.Chem import Descriptors

def load_sdf_as_rdmols(sdf_path: str) -> List[Chem.Mol]:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    for m in mols:
        try:
            AllChem.Compute2DCoords(m)
        except Exception:
            pass
    return mols

def calculate_molecular_properties(smiles: str) -> Union[Dict[str, Union[float, int]], None]:
    """
    Calculates several key molecular properties (descriptors) from a SMILES string.

    Args:
        smiles: A string representing the molecule in SMILES format (e.g., "CCO" for ethanol).

    Returns:
        A dictionary of calculated properties, or None if the SMILES is invalid.
    """
    # 1. Convert SMILES string to an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Error: Invalid SMILES string: {smiles}")
        return None

    # 3. Calculate the descriptors
    properties = {
        # General Properties
        "Molecular_Weight": Descriptors.MolWt(mol),
        
        # Lipophilicity and Polarity (important for drug absorption)
        "LogP": Descriptors.MolLogP(mol),          # Partition coefficient
        "TPSA": Descriptors.TPSA(mol),             # Topological Polar Surface Area
        
        # Structural Features
        "Num_H_Bond_Donors": Descriptors.NumHDonors(mol),
        "Num_H_Bond_Acceptors": Descriptors.NumHAcceptors(mol),

    }
    return properties
def mol_table(mols: List[Chem.Mol]) -> pd.DataFrame:
    rows = []
    for m in mols:
        name = m.GetProp("_Name") if m.HasProp("_Name") else ""
        smi = Chem.MolToSmiles(m) if m is not None else ""
        hac = m.GetNumHeavyAtoms() if m is not None else 0
        nat = m.GetNumAtoms() if m is not None else 0
        mw = Descriptors.MolWt(m) if m is not None else 0
        logp = Descriptors.MolLogP(m) if m is not None else ""
        tpsa = Descriptors.TPSA(m) if m is not None else ""
        n_H_D = Descriptors.NumHDonors(m) if m is not None else ""
        n_H_A = Descriptors.NumHAcceptors(m) if m is not None else ""
        rows.append({"name": name, 
                     "smiles": smi, 
                     "MW": mw,
                     "logP": logp,
                     "TPSA": tpsa,
                     "n_H_D": n_H_D,
                     "n_H_A": n_H_A,
                     "heavy_atoms": hac, 
                     "atoms": nat})
    return pd.DataFrame(rows)


def grid_image(mols: List[Chem.Mol], mols_per_row: int = 4, max_mols: int = 24):
    subset = mols[:max_mols]
    legends = []
    for m in subset:
        legends.append(m.GetProp("_Name") if m.HasProp("_Name") else "")
    return Draw.MolsToGridImage(subset, molsPerRow=mols_per_row, legends=legends, subImgSize=(250, 200))
