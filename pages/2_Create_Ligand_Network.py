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
from openfe_app.mapping_ui import paired_mapping_tab
from openfe_app.view3d import rdkit_mols_from_sdf
from openfe_app.network_ui import ligand_network_tab

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

st.title("Create Ligand Network")

mapping, \
network  = st.tabs(
    ["1) Paired Ligand Mapping", 
     "2) Ligand Network",]
)

with mapping:
    sdf_file = paths.inputs / "ligands.sdf"

    if not sdf_file.exists():
        st.warning("Upload ligands.sdf in Inputs first.")
        st.stop()
    else:
        mols = rdkit_mols_from_sdf(sdf_file)
        paired_mapping_tab(mols)

with network:
    sdf_file = paths.inputs / "ligands.sdf"
    if not sdf_file.exists():
        st.warning("Upload ligands.sdf in Inputs first.")
    else:
        mols = rdkit_mols_from_sdf(sdf_file)
        ligand_network_tab(mols)

st.divider()

st.markdown(
    """
This page supports **paired ligand mapping analysis and network construction**
based on structural similarity and atom mappings between ligands.

It is intended to help you:
- Define **paired ligand transformations** suitable for relative free energy (RFE) calculations
- Inspect how ligands are connected through shared chemical substructures
- Build and explore a **ligand network** that can be used for downstream simulation planning

### Workflow overview
1. **Paired Ligand Mapping**  
   Generate and inspect atom mappings between pairs of ligands.  
   These mappings define how one ligand is transformed into another and are the foundation
   for relative free energy calculations.

2. **Ligand Network**  
   Assemble ligand pairs into a network representation.  
   Nodes represent ligands, and edges represent valid transformations derived from the
   paired mappings. This view helps identify connectivity, hubs, and potential gaps in coverage.

⚠️ **Important notes**
- Ligands must be provided as an SDF file uploaded in the *Inputs* section.
- The mappings and network are based on **chemical similarity and topology**;
  they do not account for protein context or binding pose quality.
- Final network selection should be reviewed carefully before launching simulations.
"""
)
