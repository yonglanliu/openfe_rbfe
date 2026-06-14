from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_pages.select_dir import directory_picker_siderbar, directory_picker
from streamlit_pages.select_file import analysis_tsv_picker
from streamlit_pages.tab1_gather_data import tab1_design
import streamlit_pages.tab2_rbfe_network as tab2_rbfe_network
import streamlit_pages.tab3_dg_trend as tab3_dg_trend

import importlib
importlib.reload(tab2_rbfe_network)
importlib.reload(tab3_dg_trend)

# ============================================================
# Streamlit setup
# ============================================================
st.set_page_config(layout="wide")
st.title("FEP Analysis App")

# ============================================================
# Utility functions
# ============================================================
def reset_all_pickers():
    """
    Clear directory/file picker-related Streamlit session state.
    """
    for k in list(st.session_state.keys()):
        if (
            k.endswith("_current")
            or k.endswith("_show_list")
            or k.endswith("_selected")
            or "_folder_" in k
            or "picker" in k
            or "candidate" in k
        ):
            st.session_state.pop(k, None)

# ============================================================
# CSS: only color primary buttons
# Folder/file navigation buttons stay normal
# ============================================================
st.markdown(
    """
    <style>
    div.stButton > button[kind="primary"] {
        background-color: #FF8C00;
        color: white;
        border: 1px solid #FF8C00;
        border-radius: 8px;
        font-weight: 700;
    }

    div.stButton > button[kind="primary"]:hover {
        background-color: #CC7000;
        color: white;
        border: 1px solid #CC7000;
    }
    </style>
    """,
    unsafe_allow_html=True,
)




# ============================================================
# Reset button
# ============================================================
if st.button("Reset all pickers"):
    reset_all_pickers()
    st.rerun()


# ============================================================
# Tabs
# ============================================================
tab_1, tab_2, tab_3 = st.tabs(
    [
        "Gather Data",
        "RBFE Network",
        "DG Trend",
    ]
)

# ============================================================
# Select main working directory in sidebar
# ============================================================
selected_workdir = directory_picker_siderbar(
    label="Custom working directory",
    start_dir="/data/liuy48/openfe_m/simulation/",
    key_prefix="main_workdir",
)

if selected_workdir is None:
    st.info("Click **List**, navigate to a folder, then click **Use folder**.")
    st.stop()

workdir = Path(selected_workdir)


# ============================================================
# Tab 1: Gather Data
# ============================================================
with tab_1:
    tab1_design(directory_picker, workdir)

# ============================================================
# Tab 2: RBFE Network
# ============================================================
with tab_2:
    st.header("RBFE Network")

    st.info("Select an input analysis TSV file to build the RBFE network.")

    tab2_rbfe_network.tab2_design(analysis_tsv_picker, workdir)


# ============================================================
# Tab 3: DG Trend
# ============================================================
with tab_3:
    st.header("DG Trend")

    st.info("This tab can be used later for plotting dG trends.")

    with tab_3:
        tab3_dg_trend.tab3_design(
            analysis_tsv_picker=analysis_tsv_picker,
            workdir=workdir,
        )
