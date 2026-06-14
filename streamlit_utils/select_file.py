import streamlit as st
from pathlib import Path

# ============================================================
# TSV file picker
# ============================================================
def analysis_tsv_picker(
    label="Select input analysis TSV file",
    start_dir="/data/liuy48/openfe_m/simulation/",
    allowed_extensions=(".tsv",),
    key_prefix="analysis_tsv",
):
    st.subheader(label)

    current_key = f"{key_prefix}_current"
    show_key = f"{key_prefix}_show"
    candidate_key = f"{key_prefix}_candidate"
    selected_key = f"{key_prefix}_selected"

    if current_key not in st.session_state:
        st.session_state[current_key] = str(
            Path(start_dir).expanduser().resolve()
        )

    if show_key not in st.session_state:
        st.session_state[show_key] = False

    if candidate_key not in st.session_state:
        st.session_state[candidate_key] = None

    current = Path(st.session_state[current_key]).expanduser().resolve()

    st.write("Current directory:")
    st.code(str(current))

    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("List", key=f"{key_prefix}_btn_list_tsv"):
            st.session_state[show_key] = True
            st.rerun()

    with col2:
        if st.button("Hide", key=f"{key_prefix}_btn_hide_tsv"):
            st.session_state[show_key] = False
            st.rerun()

    with col3:
        if st.button("Up", key=f"{key_prefix}_btn_up_tsv"):
            parent = current.parent
            if parent.exists() and parent.is_dir():
                st.session_state[current_key] = str(parent)
                st.session_state[show_key] = True
                st.rerun()

    candidate_tsv_file = st.session_state.get(candidate_key)

    if candidate_tsv_file:
        st.write("Candidate TSV file:")
        st.code(candidate_tsv_file)

        if st.button("Use TSV file", key=f"{key_prefix}_btn_use_tsv_file", type="primary"):
            st.session_state[selected_key] = candidate_tsv_file
            st.session_state[show_key] = False
            st.success(f"Selected TSV file:\n{candidate_tsv_file}")
            st.rerun()

    selected_tsv_file = st.session_state.get(selected_key)

    if selected_tsv_file:
        st.write("Selected input analysis TSV:")
        st.code(selected_tsv_file)

    if not st.session_state[show_key]:
        return selected_tsv_file

    try:
        subdirs = sorted([p for p in current.iterdir() if p.is_dir()])
        files = sorted(
            [
                p
                for p in current.iterdir()
                if p.is_file() and p.suffix.lower() in allowed_extensions
            ]
        )
    except PermissionError:
        st.error("Permission denied.")
        return selected_tsv_file
    except FileNotFoundError:
        st.error("Directory not found.")
        return selected_tsv_file

    if subdirs:
        st.write("Folders:")
        for p in subdirs:
            if st.button(f"📁 {p.name}", key=f"{key_prefix}_folder_{str(p)}"):
                st.session_state[current_key] = str(p)
                st.session_state[show_key] = True
                st.rerun()

    if files:
        st.write("TSV files:")
        for p in files:
            if st.button(f"📄 {p.name}", key=f"{key_prefix}_file_{str(p)}"):
                st.session_state[candidate_key] = str(p)
                st.session_state[show_key] = True
                st.rerun()

    if not subdirs and not files:
        st.info("No folders or TSV files here.")

    return selected_tsv_file
