"""
Created by:Yonglan Liu
Date: 2025/12/20
"""

from __future__ import annotations

import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw


# ---------- helpers ----------
def _mol_name(m: Chem.Mol, fallback: str) -> str:
    try:
        if m.HasProp("_Name"):
            n = m.GetProp("_Name").strip()
            if n:
                return n
    except Exception:
        pass
    return fallback


def _to_openfe_component(m: Chem.Mol):
    try:
        from openfe import SmallMoleculeComponent
        return SmallMoleculeComponent.from_rdkit(m)
    except Exception:
        from openfe.setup import SmallMoleculeComponent
        return SmallMoleculeComponent.from_rdkit(m)


def _get_mapper(mapper_name: str):
    if mapper_name == "LoMapAtomMapper":
        from openfe.setup import LomapAtomMapper
        return LomapAtomMapper()
    elif mapper_name == "KartografAtomMapper":
        from openfe.setup import KartografAtomMapper
        return KartografAtomMapper()
    else:
        raise ValueError(mapper_name)


def _extract_mapping(mapping_obj):
    for attr in (
        "componentA_to_componentB",
        "atom_map",
        "mapping",
        "atom_mapping",
    ):
        if hasattr(mapping_obj, attr):
            d = getattr(mapping_obj, attr)
            if isinstance(d, dict) and d:
                return {int(k): int(v) for k, v in d.items()}
    raise RuntimeError("Cannot extract atom mapping from OpenFE mapping object.")


def _render_py3dmol_or_html(obj, height=720):
    if hasattr(obj, "_make_html"):
        st.components.v1.html(obj._make_html(), height=height)
    elif isinstance(obj, str):
        st.components.v1.html(obj, height=height)
    elif isinstance(obj, dict) and "html" in obj:
        st.components.v1.html(obj["html"], height=height)
    else:
        st.code(repr(obj))


# ---------- main UI ----------

def paired_mapping_tab(ligand_mols: list[Chem.Mol]):
    st.subheader("Paired Ligand Mapping (OpenFE)")

    if not ligand_mols:
        st.warning("No ligands loaded.")
        return

    names = [_mol_name(m, f"lig_{i}") for i, m in enumerate(ligand_mols)]

    tab_map, tab_map3d = st.tabs(["Mapping", "Mapping 3D"])

    # =======================
    # Mapping (2D + table)
    # =======================
    with tab_map:
        c1, c2, c3, c4 = st.columns([2.2, 2.2, 2.0, 1.2])

        with c1:
            ligA = st.selectbox("Ligand A", names, index=0)
        with c2:
            ligB = st.selectbox("Ligand B", names, index=min(1, len(names) - 1))
        with c3:
            mapper_name = st.selectbox(
                "Atom mapper",
                ["LoMapAtomMapper", "KartografAtomMapper"],
            )
        with c4:
            run = st.button("Compute mapping", type="primary")

        if not run:
            return

        molA = ligand_mols[names.index(ligA)]
        molB = ligand_mols[names.index(ligB)]

        compA = _to_openfe_component(molA)
        compB = _to_openfe_component(molB)
        mapper = _get_mapper(mapper_name)

        try:
            mapping = next(mapper.suggest_mappings(compA, compB))
        except StopIteration:
            st.error("No mapping suggested for this ligand pair.")
            return

        st.session_state["paired_mapping_obj"] = mapping
        st.session_state["paired_mapping_names"] = (ligA, ligB, mapper_name)

        amap = _extract_mapping(mapping)
        st.success(f"Mapped **{len(amap)} atoms** (A → B)")

        rows = []
        for a, b in sorted(amap.items()):
            rows.append(
                dict(
                    A_index=a,
                    A_elem=molA.GetAtomWithIdx(a).GetSymbol(),
                    B_index=b,
                    B_elem=molB.GetAtomWithIdx(b).GetSymbol(),
                )
            )

        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        imgA = Draw.MolToImage(molA, highlightAtoms=list(amap.keys()), size=(380, 320))
        imgB = Draw.MolToImage(molB, highlightAtoms=list(amap.values()), size=(380, 320))

        cA, cB = st.columns(2)
        with cA:
            st.image(imgA, caption=f"Ligand A: {ligA}")
        with cB:
            st.image(imgB, caption=f"Ligand B: {ligB}")

    # =======================
    # Mapping 3D (OpenFE)
    # =======================
    with tab_map3d:
        mapping = st.session_state.get("paired_mapping_obj")
        meta = st.session_state.get("paired_mapping_names")

        if not mapping:
            st.info("Compute a mapping first in the Mapping tab.")
            return

        ligA, ligB, mapper_name = meta
        st.write(f"3D mapping: **{ligA} ↔ {ligB}** using **{mapper_name}**")

        show_ids = st.toggle("Show 3D mapping", value=True)
        col_view, col_table = st.columns([3, 1], gap="medium")

        from openfe.utils import visualization_3D

        with col_view:
            view = visualization_3D.view_mapping_3d(mapping, show_atomIDs=show_ids)
            _render_py3dmol_or_html(view)
