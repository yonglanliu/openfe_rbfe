import streamlit as st
from pathlib import Path

# ============================================================
#               Working directory picker
# ============================================================
def directory_picker_siderbar(
    label,
    start_dir="./",
    key_prefix="workdir",
):
    st.sidebar.subheader(label)

    current_key = f"{key_prefix}_current"
    show_key = f"{key_prefix}_show_list"
    selected_key = f"{key_prefix}_selected"

    if current_key not in st.session_state:
        st.session_state[current_key] = str(
            Path(start_dir).expanduser().resolve()
        )

    if show_key not in st.session_state:
        st.session_state[show_key] = False

    current = Path(st.session_state[current_key]).expanduser().resolve()

    st.sidebar.code(str(current))

    col1, col2, col3, col4 = st.sidebar.columns(4)

    with col1:
        if st.button("List", key=f"{key_prefix}_btn_list"):
            st.session_state[show_key] = True
            st.rerun()

    with col2:
        if st.button("Hide", key=f"{key_prefix}_btn_hide"):
            st.session_state[show_key] = False
            st.rerun()

    with col3:
        if st.button("Up", key=f"{key_prefix}_btn_up"):
            parent = current.parent
            if parent.exists() and parent.is_dir():
                st.session_state[current_key] = str(parent)
                st.session_state[show_key] = True
                st.rerun()

    with col4:
        if st.button("Use folder", key=f"{key_prefix}_btn_use", type="primary"):
            st.session_state[selected_key] = str(current)
            st.session_state[show_key] = False
            st.sidebar.success(f"Selected:\n{current}")
            st.rerun()

    selected_dir = st.session_state.get(selected_key)


    if not st.session_state[show_key]:
        return selected_dir

    try:
        subdirs = sorted([p for p in current.iterdir() if p.is_dir()])
    except PermissionError:
        st.sidebar.error("Permission denied.")
        return selected_dir
    except FileNotFoundError:
        st.sidebar.error("Directory not found.")
        return selected_dir

    if not subdirs:
        st.sidebar.info("No subdirectories here.")
        return selected_dir

    st.sidebar.write("Folders:")

    for p in subdirs:
        if st.sidebar.button(
            f"📁 {p.name}",
            key=f"{key_prefix}_folder_{str(p)}",
        ):
            st.session_state[current_key] = str(p)
            st.session_state[show_key] = True
            st.rerun()

    return selected_dir


def directory_picker(
    label,
    start_dir="./",
    key_prefix="workdir",
):
    st.subheader(label)

    current_key = f"{key_prefix}_current"
    show_key = f"{key_prefix}_show_list"
    selected_key = f"{key_prefix}_selected"

    if current_key not in st.session_state:
        st.session_state[current_key] = str(
            Path(start_dir).expanduser().resolve()
        )

    if show_key not in st.session_state:
        st.session_state[show_key] = False

    current = Path(st.session_state[current_key]).expanduser().resolve()

    st.code(str(current))
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        if st.button("List", key=f"{key_prefix}_btn_list"):
            st.session_state[show_key] = True
            st.rerun()

    with col2:
        if st.button("Hide", key=f"{key_prefix}_btn_hide"):
            st.session_state[show_key] = False
            st.rerun()

    with col3:
        if st.button("Up", key=f"{key_prefix}_btn_up"):
            parent = current.parent
            if parent.exists() and parent.is_dir():
                st.session_state[current_key] = str(parent)
                st.session_state[show_key] = True
                st.rerun()

    with col4:
        if st.button("Use folder", key=f"{key_prefix}_btn_use", type="primary"):
            st.session_state[selected_key] = str(current)
            st.session_state[show_key] = False
            st.success(f"Selected:\n{current}")
            st.rerun()

    selected_dir = st.session_state.get(selected_key)

    # if selected_dir:
    #     st.write("Selected directory:")
    #     st.code(selected_dir)

    if not st.session_state[show_key]:
        return selected_dir

    try:
        subdirs = sorted([p for p in current.iterdir() if p.is_dir()])
    except PermissionError:
        st.error("Permission denied.")
        return selected_dir
    except FileNotFoundError:
        st.error("Directory not found.")
        return selected_dir

    if not subdirs:
        st.info("No subdirectories here.")
        return selected_dir

    st.write("Folders:")

    for p in subdirs:
        if st.button(
            f"📁 {p.name}",
            key=f"{key_prefix}_folder_{str(p)}",
        ):
            st.session_state[current_key] = str(p)
            st.session_state[show_key] = True
            st.rerun()

    return selected_dir
