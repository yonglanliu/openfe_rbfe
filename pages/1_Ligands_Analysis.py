"""
Created by: Yonglan Liu
Date: 2025/12/27
"""

from __future__ import annotations

import streamlit as st
from pathlib import Path
import pandas as pd
import os

from openfe_app.config import ProjectPaths
from openfe_app.io_utils import save_upload
from rdkit import Chem
from openfe_app.ligands_viz import (
    load_sdf_as_rdmols, 
    grid_image, 
    mol_table
)
from openfe_app.ligands_viz import load_sdf_as_rdmols, grid_image, mol_table
from openfe_app.view3d import (
    extract_ligands_from_pdb,
    rdkit_mol_from_pdb_block,
    rdkit_mols_from_sdf,
    align_to_reference,
    mol_name,
    render_protein_ligands_single_view, 
)
if "project_root" not in st.session_state:
    st.warning("No project selected yet. Go to the Project page and create/select one.")
    st.stop()
paths = ProjectPaths(Path(st.session_state["project_root"]))
st.session_state["project_root"] = str(paths.root)

st.sidebar.write("Folders:")
st.sidebar.code(
    f"inputs:   {paths.inputs}\n"
    f"prepared:  {paths.planned}\n"
    f"results:  {paths.results}\n"
    f"analysis: {paths.analysis}"
)
st.title("Ligand Analysis")

tab_inputs, \
tab_ligand_structures, \
tab_3d = st.tabs(
    ["1) Upload", 
     "2) Ligand Structures", 
     "3) Binding Mode"]
)

# ---------- 1) upload ----------
with tab_inputs:
    st.subheader("Upload inputs")
    st.write("Upload an SDF of ligands and a prepared protein PDB.")

    col1, col2 = st.columns(2)

    with col1:
        sdf_up = st.file_uploader("Ligands SDF (.sdf)", type=["sdf"])
        if sdf_up:
            sdf_path = save_upload(sdf_up, paths.inputs / "ligands.sdf")
            st.success(f"Saved: {sdf_path}")

    with col2:
        pdb_up = st.file_uploader("Protein PDB (.pdb)", type=["pdb"])
        if pdb_up:
            pdb_path = save_upload(pdb_up, paths.inputs / "protein.pdb")
            st.success(f"Saved: {pdb_path}")
    
    st.divider()

    st.markdown(
    """
    This app supports **ligand structure visualization and binding mode analysis**.

    You can upload ligand SDF files and a prepared protein PDB, align ligands to a reference,
    and inspect binding-site interactions in an interactive 3D viewer.


    ⚠️ **Important note:**  
    This tool performs **structural alignment only**. Ligands are *not docked* into the protein;
    their placement reflects alignment to an existing reference ligand or to each other.
    """
    )

# ---------- 2) Ligands ----------
with tab_ligand_structures:
    st.subheader("Ligand viewer (RDKit)")

    sdf_file = paths.inputs / "ligands.sdf"
    if not sdf_file.exists():
        st.warning("Upload ligands.sdf in the Inputs tab first.")
    else:
        mols = load_sdf_as_rdmols(str(sdf_file))
        st.write(f"Loaded **{len(mols)}** ligands.")


        max_mols = st.slider("Max ligands to display", 4, 60, 24, step=4)
        per_row = st.slider("Molecules per row", 2, 6, 4)

        img = grid_image(mols, mols_per_row=per_row, max_mols=max_mols)
        st.image(img)

        st.divider()
        st.subheader("Molecule Information")
        df = mol_table(mols)
        st.dataframe(df, use_container_width=True)

# ---------- 3) view 3D structure ----------
with tab_3d:
    st.subheader("3D: Align ligands and view with protein")

    sdf_file = paths.inputs / "ligands.sdf"
    pdb_file = paths.inputs / "protein.pdb"

    if not sdf_file.exists():
        st.warning("Upload ligands.sdf in Inputs first.")
        st.stop()
    if not pdb_file.exists():
        st.warning("Upload protein.pdb in Inputs first.")
        st.stop()

    # Extract ligands from sdf file
    mols = rdkit_mols_from_sdf(sdf_file)
    ligand_names = [mol_name(ref) for ref in mols]
    if len(mols) == 0:
        st.error("No valid molecules found in ligands.sdf")
        st.stop()

    st.info(f"Loaded **{len(mols)}** ligands from SDF")

    # --- Find a reference ligand pose from the protein PDB (if present) ---
    ligs_in_pdb = extract_ligands_from_pdb(pdb_file)
    ref_ligand_rdkit = None
    ref_label = None

    if ligs_in_pdb:
        st.info(
            f"Found {len(ligs_in_pdb)} non-protein residues in PDB. "
            f"Will use the largest as reference pose: {ligs_in_pdb[0].resname} "
            f"(chain {ligs_in_pdb[0].chain}, resid {ligs_in_pdb[0].resid})."
        )
        ref_ligand_rdkit = rdkit_mol_from_pdb_block(ligs_in_pdb[0].pdb_block)
        ref_label = f"PDB ligand {ligs_in_pdb[0].resname}:{ligs_in_pdb[0].chain}{ligs_in_pdb[0].resid}"
        if ref_ligand_rdkit is None:
            st.warning(
                "Could not parse PDB ligand into RDKit (common if bond orders are ambiguous). "
                "Will fallback to ligand-to-ligand alignment only."
            )
    else:
        st.warning(
            "No ligand detected in the protein PDB (no suitable HETATM group). "
            "I can still align ligands to each other, but placement in the binding site is NOT known "
            "(this is not docking)."
        )

    # --- UI controls ---
    col1, col2, col3, col4 = st.columns([2, 2, 2, 2])
    with col1:
        max_show = st.slider("Max ligands to show in 3D", 1, min(100, len(mols)), min(20, len(mols)))
    with col2:
        align_mode = st.selectbox(
            "Alignment reference",
            options=ligand_names + ([ref_label] if ref_ligand_rdkit is not None else []),
            # options=["First ligand in SDF"] + ([ref_label] if ref_ligand_rdkit is not None else []),
        )
    with col3:
        show_style = st.selectbox("Ligand style", ["stick", "sphere", "line"], index=0)
    with col4:
        show_labels = st.toggle("Show panel labels", value=False)

    subset = mols[:max_show]
    #st.text(subset)

    # Choose reference
    if align_mode in ligand_names:
        index_of_ligand = ligand_names.index(align_mode)
        ref = subset[index_of_ligand]
    else:
        ref = ref_ligand_rdkit  # type: ignore

    aligned = align_to_reference(subset, ref)

    # --- Read protein PDB text ---
    with open(pdb_file, "r") as f:
        protein_pdb = f.read()

    # --- Render ONE window with options (protein + ALL aligned ligands) ---
    aligned_sdf_blocks = [Chem.MolToMolBlock(m) for m in aligned]
    aligned_names = [mol_name(m, f"lig_{i}") for i, m in enumerate(aligned)]

    render_protein_ligands_single_view(
        protein_pdb_text=protein_pdb,
        aligned_ligands=aligned_sdf_blocks,   # list[str] MolBlocks
        ligand_names=aligned_names,
        key_prefix="tab3d_single",
        width=1100,
        height=680,
    )

    st.caption(
        "Note: This view overlays ligands by 3D alignment. If the PDB contains a bound ligand, "
        "ligands are aligned to that pose; otherwise they are aligned to the ligand you selected and shown near the protein "
        "(not a docking result)."
    )

st.divider()


