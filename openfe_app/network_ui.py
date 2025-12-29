from __future__ import annotations

import streamlit as st
from rdkit import Chem
import pandas as pd
import io

def _ensure_component_names(ligands):
    """Ensure each OpenFE component has a stable name (used in dropdowns)."""
    # Most OpenFE components expose .name; we set it from rdkit _Name if needed.
    for i, comp in enumerate(ligands):
        nm = getattr(comp, "name", None)
        if nm is None or str(nm).strip() == "":
            # fallback
            try:
                comp.name = f"lig_{i}"
            except Exception:
                pass


def _to_components(rdmols):
    import openfe

    comps = []
    for i, m in enumerate(rdmols):
        if not m.HasProp("_Name") or not m.GetProp("_Name").strip():
            m.SetProp("_Name", f"lig_{i}")
        comps.append(openfe.SmallMoleculeComponent.from_rdkit(m))
    _ensure_component_names(comps)
    return comps


def _get_mapper(name: str):
    from openfe import setup
    if name == "LomapAtomMapper":
        return setup.LomapAtomMapper()
    if name == "KartografAtomMapper":
        return setup.KartografAtomMapper()
    raise ValueError(name)


def _get_planner(name: str):
    import openfe
    planners = openfe.ligand_network_planning

    if name == "generate_minimal_spanning_network":
        return planners.generate_minimal_spanning_network, False
    if name == "generate_minimal_redundant_network":
        return planners.generate_minimal_redundant_network, False
    if name == "generate_radial_network":
        return planners.generate_radial_network, True

    raise ValueError(name)


def ligand_network_tab(rdmols):
    st.subheader("Ligand Network (OpenFE)")

    if not rdmols:
        st.warning("No ligands loaded.")
        return

    # Controls
    c1, c2, c3 = st.columns([2, 2, 1.2])
    with c1:
        mapper_name = st.selectbox("Atom mapper", ["LomapAtomMapper", "KartografAtomMapper"], index=0)
    with c2:
        planner_name = st.selectbox(
            "Network planner",
            [
                "generate_minimal_spanning_network",
                "generate_minimal_redundant_network",
                "generate_radial_network",
            ],
            index=0,
        )
    with c3:
        build = st.button("Build network", type="primary")

    # Build network only when requested (avoid slow rebuild on every widget change)
    if build:
        import openfe
        from openfe import lomap_scorers

        ligands = _to_components(rdmols)
        mapper = _get_mapper(mapper_name)
        planner_fn, needs_central = _get_planner(planner_name)

        scorer = lomap_scorers.default_lomap_score

        central_ligand = None
        if needs_central:
            # show selector ONLY for radial
            central_name = st.selectbox(
                "Central ligand (radial network)",
                options=[c.name for c in ligands],
                index=0,
                key="net_radial_center",
            )
            name_to_comp = {c.name: c for c in ligands}
            central_ligand = name_to_comp[central_name]

            ligand_network = planner_fn(
                ligands=ligands,
                central_ligand=central_ligand,
                mappers=[mapper],
                scorer=scorer,
            )
        else:
            ligand_network = planner_fn(
                ligands=ligands,
                mappers=[mapper],
                scorer=scorer,
            )

        st.session_state["ligand_network_obj"] = ligand_network
        st.session_state["ligand_network_meta"] = {
            "mapper": mapper_name,
            "planner": planner_name,
            "central": getattr(central_ligand, "name", None) if central_ligand else None,
            "n_nodes": len(ligand_network.nodes),
            "n_edges": len(ligand_network.edges),
        }

        meta = st.session_state["ligand_network_meta"]
        extra = f", central={meta['central']}" if meta["central"] else ""
        st.success(f"Built network: **{meta['n_nodes']} nodes**, **{meta['n_edges']} edges**{extra}")

    ligand_network = st.session_state.get("ligand_network_obj")
    meta = st.session_state.get("ligand_network_meta")

    if not ligand_network:
        st.info("Click **Build network** to generate a LigandNetwork.")
        return

    # --- Plot controls (make it smaller) ---
    st.markdown("### Network")
    pc1, pc2, pc3 = st.columns([1.2, 1.2, 2])
    with pc1:
        fig_w = st.slider("Plot width (in)", 4.0, 14.0, 7.0, 0.5, key="net_fig_w")
    with pc2:
        fig_h = st.slider("Plot height (in)", 3.0, 12.0, 5.0, 0.5, key="net_fig_h")
    with pc3:
        st.write(f"Planner: **{meta['planner']}**, Mapper: **{meta['mapper']}**" + (f", Central: **{meta['central']}**" if meta.get("central") else ""))

    # --- Visualize topology ---
    from openfe.utils.atommapping_network_plotting import plot_atommapping_network

    fig = plot_atommapping_network(ligand_network)
    try:
        fig.set_size_inches(fig_w, fig_h)
        fig.tight_layout()
    except Exception:
        pass
    col_view, col_table = st.columns([1, 1], gap="medium")

    with col_view:
        st.pyplot(fig, use_container_width=False)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
        buf.seek(0)

        st.download_button(
            label="Download figure (PNG)",
            data=buf,
            file_name="figure.png",
            mime="image/png",
        )

    with col_table:
        edges = list(ligand_network.edges)
        st.markdown("### Edges")
        st.caption(f"Total edges: {len(edges)}")

        if edges:
            df_edges = pd.DataFrame(
                [
                    {
                        "Ligand A": e.componentA.name,
                        "Ligand B": e.componentB.name,
                    }
                    for e in edges
                ]
            )

            st.dataframe(
                df_edges.head(200),
                use_container_width=True,
                height=600,
                hide_index=True,
            )

            if len(edges) > 200:
                st.info("Showing first 200 edges.")
        else:
            st.info("No edges found.")

