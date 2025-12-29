from __future__ import annotations

import json
import streamlit as st
from pathlib import Path
import pandas as pd
import os
from io import BytesIO
from openfe_app.config import ProjectPaths
from openfe_app.io_utils import save_upload
from openfe_app.json_analysis import (
    extract_dg,
    extract_f_and_r_convergence, 
    extract_rex,
    extract_mbar
)
from openfe_app.streamlit_plot import (
    plot_fr_convergence,
    plot_replica_exchange_matrix,
    plot_mbar_matrix,
    fig_to_png_bytes
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
replicates = os.listdir(paths.results)
replicate_sel = st.sidebar.selectbox("Select Replicate", replicates, index=0)
st.session_state["replcate_sel"] = replicate_sel
replicate_path = paths.results / replicate_sel
output_jsons = replicate_path.glob("*.json")
if output_jsons is not None:
    base_names = set()
    output_json_files = []
    for name in output_jsons:
        bn = os.path.basename(name)
        output_json_files.append(str(bn))
        bn = bn.replace(".json", "")
        if "solvent" in bn:
            bn = bn.replace("_solvent_", "_to_")
            bn = bn.replace("_solvent", "")
        if "complex" in bn:
            bn = bn.replace("_complex_", "_to_")
            bn = bn.replace("_complex", "")
        base_names.add(bn)
    st.session_state["output_json_files"] = output_json_files
else:
    st.sidebar.info(f"There isn't JSON file in folder {replicate_path}")

sel_transformation = st.sidebar.selectbox("Select Transformation", list(base_names), index=0)
st.session_state['selected_transformation'] = sel_transformation

nameA = st.session_state['selected_transformation'].split("_to_")[0]
nameB = st.session_state['selected_transformation'].split("_to_")[1]

sol_name = nameA + "_solvent_" + nameB + "_solvent" + ".json"
complex_name = nameA + "_complex_" + nameB + "_complex" + ".json"
sol_json_file = Path(paths.results / st.session_state["replcate_sel"]/ sol_name)
complex_json_file = Path(paths.results / st.session_state["replcate_sel"]/ complex_name)

# -----------------------------------------------------
st.set_page_config(page_title="OpenFE Dashboard", layout="wide")
st.title("Analysis for single transformation")
# ---------- Tabs ----------
tab_intro, final_dg, f_r_convergence, rex, mbar = st.tabs(
    ["1) Introduction", "2) dG data", "3) forward/reverse convergence", "4) replica exchange", "5) MBAR analysis"]
)

# ------------------------------------------------
#                 Introduction
# ------------------------------------------------
with tab_intro:
    st.markdown(f"""
# 📊 Results Analysis (RBFE)

This page helps you **inspect, validate, and visualize** the outputs of an OpenFE **Relative Binding Free Energy (RBFE)** calculation.
It reads the OpenFE JSON result files for a chosen transformation and summarizes key diagnostics for both legs:

- **Solvent leg** (`{paths.results}/replicate_[index]/*_solvent.json`)
- **Complex leg** (`{paths.results}/replicate_[index]/*_complex.json`)

---

## ✅ How to use this page

1. **Select a replicate and a transformation** in the left sidebar (e.g., `replicate_0` and `ligA_to_ligB`).
2. The dashboard will automatically load the corresponding JSON results from:
   - `{paths.results}/replicate_[index]/*_solvent.json`
   - `{paths.results}/replicate_[index]/*_complex.json`
3. Navigate through the tabs above to explore different analyses and download figures.

---

## 🧭 What each tab shows

### 2) ΔG data
Displays per-repeat estimates of **ΔG** for each leg (solvent / complex).  
Use this to quickly spot outliers or inconsistencies across repeats.

### 3) Forward/Reverse convergence
Plots **forward and reverse cumulative estimates** to help judge convergence stability.  
Good runs typically show forward and reverse curves approaching a stable plateau.

### 4) Replica exchange
Shows the **replica exchange transition matrix** across λ-windows.  
A healthy run typically has transitions spread around the diagonal (not stuck).

### 5) MBAR analysis
Visualizes the **MBAR overlap / weight matrix**, a useful diagnostic for:
- window overlap quality  
- whether λ spacing is adequate  
- identifying poorly connected regions

---

## ⚠️ Common interpretation tips

- If **ΔG varies strongly across repeats**, you may need longer sampling or better λ spacing.
- If replica exchange is **mostly zeros** or very sparse, replicas may be stuck (bad mixing).
- If MBAR overlap looks **blocky / disconnected**, adjacent windows may not overlap enough.

---

## 🔗 References

- OpenFE protocol configuration:  
  [Choose and Configure a Protocol](https://docs.openfree.energy/en/stable/cookbook/choose_protocol.html)

- OpenFE documentation home:  
  [OpenFE Docs](https://docs.openfree.energy/)
""")

    st.info("Tip: Start with **ΔG data**, then check **Forward/Reverse convergence** before deeper diagnostics.")

st.divider()

# --------------------------------------------
with open(sol_json_file, "r") as f:
    st.session_state["results_solvent_json"] = json.load(f)
with open(complex_json_file, "r") as f:
    st.session_state["results_complex_json"] = json.load(f)

#st.text(st.session_state["results_solvent_json"])
# --------------------------------------------------
with final_dg:
    st.markdown(f"## $\Delta G$ data")
    col1, col2 = st.columns([2, 2])

    # solvation leg
    with col1:
        st.markdown("### Solvent Leg")
        data_sol = st.session_state["results_solvent_json"]
        df_sol = extract_dg(data_sol)
        st.dataframe(df_sol, use_container_width=True)
    with col2:
        st.markdown("### Complex Leg")
        data_complex = st.session_state["results_complex_json"]
        df_complex = extract_dg(data_complex)
        st.dataframe(df_complex, use_container_width=True)
    st.divider()
    c1, c2, c3 = st.columns([2, 2, 2])
    with c2:
        ddG = round(df_complex["DG (kcal/mol)"].values[-1] - df_sol["DG (kcal/mol)"].values[-1], 2)
        sigma = (
        df_complex["uncertainty (kcal/mol)"].values[-1] ** 2
        + df_sol["uncertainty (kcal/mol)"].values[-1] ** 2
        ) ** 0.5
        st.markdown(
                    r"""
                    $$
                    \Delta\Delta G_{\mathrm{bind}}
                    =
                    \Delta G_{\mathrm{complex}} - \Delta G_{\mathrm{solvent}}
                    $$
                    """
                    )


        st.latex(rf"\Delta\Delta G = {ddG:.2f}\ \mathrm{{kcal/mol}}")
        st.divider()
        st.markdown(
                    r"""
                    $$
                    \sigma_{\Delta\Delta G}
                    =
                    \sqrt{\sigma_{\mathrm{complex}}^2 + \sigma_{\mathrm{solvent}}^2}
                    $$
                    """
                    )

        st.latex(rf"\sigma_{{\Delta\Delta G}} = {sigma:.2f}\ \mathrm{{kcal/mol}}")


with f_r_convergence:
    st.subheader("Forward/Reverse Convergency")
    st.divider()
    colA, colB, colC = st.columns(3)
    with colA:
        show_errorbars = st.checkbox("Show error bars", value=True)
        marker_size = st.slider("Marker size", 3, 14, 6)
        line_width = st.slider("Line width", 1.0, 5.0, 2.5, 0.1)
    with colB:
        show_band = st.checkbox("Show horizontal band", value=True)
        band_halfwidth = st.slider("Band halfwidth (kcal/mol)", 0.0, 1.0, 0.10, 0.01)
        band_alpha = st.slider("Band alpha", 0.0, 1.0, 0.15, 0.01)
    with colC:
        band_mode = st.radio("Band center", ["Use last forward point", "Use final estimate"], horizontal=False)

    st.divider()
    data_sol = st.session_state["results_solvent_json"]
    data_complex = st.session_state["results_complex_json"]
    
    repeats_sol = list(data_sol["unit_results"].keys())
    repeats_com = list(data_complex["unit_results"].keys())
    assert len(repeats_sol) == len(repeats_com)
    for i in range(len(repeats_sol)):
        rep_s = data_sol["unit_results"][repeats_sol[i]]
        rep_c = data_complex["unit_results"][repeats_com[i]]

        col1,_, col2 = st.columns([2,1,2])

        with col1:
            st.write(f"Solvent_repeat_{i}")
            df_s = extract_f_and_r_convergence(rep_s)
            fig = plot_fr_convergence(df_s,
                                      marker_size=marker_size,
                                      line_width=line_width,
                                      band_halfwidth=band_halfwidth,
                                      alpha_band=band_alpha)
            png_bytes = fig_to_png_bytes(fig, dpi=300)

            st.pyplot(fig, use_container_width=True)
            st.download_button(
                    label="Download figure (PNG)",
                    data=png_bytes,
                    file_name=f"solvent_forward_reverse_repeat{i}.png",
                    mime="image/png",)
            st.dataframe(df_s, use_container_width=True)
        with col2:
            st.write(f"Complex_repeat_{i}")
            df_c = extract_f_and_r_convergence(rep_c)
            fig = plot_fr_convergence(df_c, 
                                      marker_size=marker_size,
                                      line_width=line_width,
                                      band_halfwidth=band_halfwidth,
                                      alpha_band=band_alpha)
            st.pyplot(fig, use_container_width=True)
            png_bytes = fig_to_png_bytes(fig, dpi=300)
            st.download_button(
                label="Download figure (PNG)",
                data=png_bytes,
                file_name=f"complex_forward_reverse_repeat{i}.png",
                mime="image/png",
            )
            st.dataframe(df_c, use_container_width=True)

        st.divider()   # visual row separator

with rex:
    st.subheader("Replica Exchange Transaction")
    st.divider()
    # ----------------------------------------------------------
    data_sol = st.session_state["results_solvent_json"]
    data_complex = st.session_state["results_complex_json"]
    # ----------------------------------------------------------

    repeats_sol = list(data_sol["unit_results"].keys())
    repeats_com = list(data_complex["unit_results"].keys())
    assert len(repeats_sol) == len(repeats_com)
    for i in range(len(repeats_sol)):
        rep_s = data_sol["unit_results"][repeats_sol[i]]
        rep_c = data_complex["unit_results"][repeats_com[i]]

        col1,_, col2 = st.columns([2,1,2])

        with col1:
            st.write(f"Solvent_repeat_{i}")
            m, _ = extract_rex(rep_s)
            m_plot = m.copy()
            m_plot[m_plot<0.01] = 0
            fig = plot_replica_exchange_matrix(m_plot)
            png_bytes = fig_to_png_bytes(fig, dpi=300)

            st.pyplot(fig, use_container_width=True)
            st.download_button(
                    label="Download figure (PNG)",
                    data=png_bytes,
                    file_name=f"solvent_forward_reverse_repeat{i}.png",
                    mime="image/png",
                    key=f"download_solvent_rex_{i}")
        with col2:
            st.write(f"Complex_repeat_{i}")
            m, _ = extract_rex(rep_c)
            m_plot = m.copy()
            m_plot[m_plot<0.01] = 0
            fig = plot_replica_exchange_matrix(m_plot)
            png_bytes = fig_to_png_bytes(fig, dpi=300)
            st.pyplot(fig, use_container_width=True)
            st.download_button(
                label="Download figure (PNG)",
                data=png_bytes,
                file_name=f"replica_exchange_transaction_repeat{i}.png",
                mime="image/png",
                key=f'download_complex_rex_{i}'
            )

        st.divider()   # visual row separator

with mbar:
    st.subheader("Mbar")
    st.divider()
    # ----------------------------------------------------------
    data_sol = st.session_state["results_solvent_json"]
    data_complex = st.session_state["results_complex_json"]
    # ----------------------------------------------------------

    repeats_sol = list(data_sol["unit_results"].keys())
    repeats_com = list(data_complex["unit_results"].keys())
    assert len(repeats_sol) == len(repeats_com)
    for i in range(len(repeats_sol)):
        rep_s = data_sol["unit_results"][repeats_sol[i]]
        rep_c = data_complex["unit_results"][repeats_com[i]]

        col1, _, col2 = st.columns([2,1,2])

        with col1:
            st.write(f"Solvent_repeat_{i}")
            m, _ = extract_mbar(rep_s)
            m_plot = m.copy()
            m_plot[m_plot<0.01] = 0
            fig = plot_mbar_matrix(m_plot)
            png_bytes = fig_to_png_bytes(fig, dpi=300)

            st.pyplot(fig, use_container_width=True)
            st.download_button(
                    label="Download figure (PNG)",
                    data=png_bytes,
                    file_name=f"mbar_repeat{i}.png",
                    mime="image/png",
                    key=f"download_solvent_mbar_{i}")
        with col2:
            st.write(f"Complex_repeat_{i}")
            m, _ = extract_mbar(rep_c)
            m_plot = m.copy()
            m_plot[m_plot<0.01] = 0
            fig = plot_mbar_matrix(m_plot)
            png_bytes = fig_to_png_bytes(fig, dpi=300)
            st.pyplot(fig, use_container_width=True)
            st.download_button(
                label="Download figure (PNG)",
                data=png_bytes,
                file_name=f"mbar_repeat{i}.png",
                mime="image/png",
                key=f'download_complex_mbar_{i}'
            )

        st.divider()   # visual row separator