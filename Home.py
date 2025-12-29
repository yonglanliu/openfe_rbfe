from __future__ import annotations

import streamlit as st
from pathlib import Path
import os

from openfe_app.config import ProjectPaths

st.set_page_config(page_title="OpenFE Dashboard", layout="wide")

st.title("Alchemical FEP Calculation Using OpenFE")
st.header("FEP Project")

project_name = st.text_input("Project Name")
# If user changes project name, reset the saved project_root
prev_name = st.session_state.get("prev_project_name")

if project_name and prev_name != project_name:
    st.session_state.pop("project_root", None)

st.session_state["prev_project_name"] = project_name
default_root = str(Path.cwd() / project_name) if project_name else str(Path.cwd())

if st.button("Create Project"):
    project_root = st.sidebar.text_input(
    "Project directory",
    value=default_root,
    )
    paths = ProjectPaths(Path(project_root).resolve())
    paths.ensure()
    st.session_state["project_root"] = str(paths.root)

    st.sidebar.write("Folders:")
    st.sidebar.code(
        f"inputs:   {paths.inputs}\n"
        f"planned:  {paths.planned}\n"
        f"results:  {paths.results}\n"
        f"analysis: {paths.analysis}"
    )
    st.session_state["project_root"] = str(paths.root)


st.divider()

st.markdown("""# OpenFE RBFE Dashboard
An intuitive, web-based dashboard for performing Relative Binding Free Energy (RBFE) and alchemical Free Energy Perturbation (FEP) calculations using the [OpenFE framework]("https://github.com/IIIS-Li-Group/OpenFE").

## 📖 Background & Concepts
If you are new to the world of alchemical transformations, start with my introductory guide: 
👉 Read the Blog Post: [Introduction to Alchemical Free Energy Calculation: FEP and TI for RBFE](https://yonglanliu.github.io/2025/12/19/RBFE.html). 

This guide breaks down the core concepts without overwhelming you with complex mathematical formulas.

---

## 🚀 Key Features
This package streamlines the standard OpenFE workflow into five user-friendly steps:
1. Molecular Visualization: Interactively view your protein target and ligand structures within the dashboard.
2. Network Planning: Automatically generate and inspect the Ligand Network (transformation edges).
3. System Setup: Configure your alchemical system, including solvation, force fields, and box parameters.
4. Job Execution: Simplified controls to prepare and run simulation scripts.
5. Results Analysis: Parse raw output data, visualize convergence, and extract final $\Delta \Delta G$ values.

---
""")



