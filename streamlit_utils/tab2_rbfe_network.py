import os
import json
import time
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import streamlit as st
import streamlit.components.v1 as components
from pyvis.network import Network


# ============================================================
# Default column names
# ============================================================
DEFAULT_DG_COL = "DG(i->j) (kcal/mol)"
DEFAULT_UNC_COL = "MBAR uncertainty (kcal/mol)"


# ============================================================
# Helper functions
# ============================================================

def add_persist_node_positions_js(html, storage_key="rbfe_node_positions"):
    """
    Save dragged PyVis node positions in browser localStorage and restore them
    after Streamlit reruns/rebuilds the graph.
    """

    persist_js = f"""
    <script>
    function restoreNodePositions() {{
        try {{
            const saved = localStorage.getItem("{storage_key}");
            if (!saved) return;

            const positions = JSON.parse(saved);

            if (typeof nodes === "undefined" || typeof network === "undefined") {{
                setTimeout(restoreNodePositions, 200);
                return;
            }}

            Object.keys(positions).forEach(function(nodeId) {{
                nodes.update({{
                    id: nodeId,
                    x: positions[nodeId].x,
                    y: positions[nodeId].y,
                    physics: false,
                    fixed: false
                }});
            }});

            network.redraw();
        }} catch (err) {{
            console.log("Could not restore node positions:", err);
        }}
    }}

    function saveNodePositions() {{
        try {{
            if (typeof network === "undefined") return;

            const positions = network.getPositions();
            localStorage.setItem("{storage_key}", JSON.stringify(positions));
            console.log("Saved node positions:", positions);
        }} catch (err) {{
            console.log("Could not save node positions:", err);
        }}
    }}

    function attachPositionSaver() {{
        if (typeof network === "undefined") {{
            setTimeout(attachPositionSaver, 200);
            return;
        }}

        restoreNodePositions();

        network.on("dragEnd", function(params) {{
            saveNodePositions();
        }});

        network.on("stabilized", function() {{
            restoreNodePositions();
        }});
    }}

    setTimeout(attachPositionSaver, 500);
    </script>
    """

    return html.replace("</body>", persist_js + "\n</body>")

def weighted_mean_and_uncertainty(values, uncertainties):
    """
    Inverse-variance weighted mean and uncertainty.
    """
    v = np.asarray(values, dtype=float)
    s = np.asarray(uncertainties, dtype=float)

    if np.any(s <= 0) or np.any(~np.isfinite(s)):
        raise ValueError("All uncertainties must be positive and finite.")

    w = 1.0 / s**2
    mean = np.sum(w * v) / np.sum(w)
    unc = np.sqrt(1.0 / np.sum(w))

    return mean, unc


def build_fep_graph(df, dg_col, unc_col):
    """
    Build directed FEP graph from replicate-level FEP results.
    """
    rows = []

    for (leg, lig_i, lig_j), g in df.groupby(["leg", "ligand_i", "ligand_j"]):
        dg_mean, dg_unc = weighted_mean_and_uncertainty(
            g[dg_col],
            g[unc_col],
        )

        rows.append(
            {
                "leg": str(leg),
                "ligand_i": str(lig_i),
                "ligand_j": str(lig_j),
                "DG_mean": dg_mean,
                "DG_unc": dg_unc,
                "n_replicates": len(g),
            }
        )

    summary = pd.DataFrame(rows)

    G = nx.DiGraph()

    ligands = pd.unique(df[["ligand_i", "ligand_j"]].values.ravel("K"))

    for lig in ligands:
        G.add_node(str(lig))

    for (lig_i, lig_j), sub in summary.groupby(["ligand_i", "ligand_j"]):
        lig_i = str(lig_i)
        lig_j = str(lig_j)

        attrs = {
            "ligand_i": lig_i,
            "ligand_j": lig_j,
        }

        for _, row in sub.iterrows():
            leg = row["leg"]

            attrs[f"{leg}_DG"] = float(row["DG_mean"])
            attrs[f"{leg}_unc"] = float(row["DG_unc"])
            attrs[f"{leg}_n"] = int(row["n_replicates"])

        if "complex_DG" in attrs and "solvent_DG" in attrs:
            attrs["ddG"] = attrs["complex_DG"] - attrs["solvent_DG"]
            attrs["ddG_unc"] = np.sqrt(
                attrs["complex_unc"] ** 2 + attrs["solvent_unc"] ** 2
            )

        G.add_edge(lig_i, lig_j, **attrs)

    return G, summary


def add_download_png_button(html):
    """
    Add a browser-side button to download the currently visible PyVis canvas as PNG.
    """
    download_js = """
    <div style="position: fixed; top: 10px; right: 10px; z-index: 9999;">
        <button onclick="downloadNetworkPNG()"
                style="font-size: 16px; padding: 8px 14px; background: white; border: 1px solid black; cursor: pointer;">
            Download PNG
        </button>
    </div>

    <script>
    function downloadNetworkPNG() {
        var canvas = document.querySelector("canvas");
        if (!canvas) {
            alert("No canvas found.");
            return;
        }

        var link = document.createElement("a");
        link.download = "fep_network_adjusted.png";
        link.href = canvas.toDataURL("image/png");
        link.click();
    }
    </script>
    """

    return html.replace("</body>", download_js + "\n</body>")


def make_initial_positions(G, layout_name, scale, rotation_deg=0):
    """
    Generate initial node positions for PyVis.

    The rotation only affects the starting layout.
    After rendering, nodes can still be dragged manually.
    """
    if layout_name == "spring":
        pos = nx.spring_layout(G, seed=42, k=1.2)
    elif layout_name == "circular":
        pos = nx.circular_layout(G)
    elif layout_name == "shell":
        pos = nx.shell_layout(G)
    elif layout_name == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42, k=1.2)

    theta = np.deg2rad(rotation_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    pyvis_pos = {}

    for node, xy in pos.items():
        x0 = float(xy[0])
        y0 = float(xy[1])

        x_rot = x0 * cos_t - y0 * sin_t
        y_rot = x0 * sin_t + y0 * cos_t

        x = x_rot * scale
        y = -y_rot * scale

        pyvis_pos[node] = {
            "x": x,
            "y": y,
        }

    return pyvis_pos


def make_edge_table(G):
    rows = []

    for u, v, d in G.edges(data=True):
        row = {
            "ligand_i": u,
            "ligand_j": v,
        }
        row.update(d)
        rows.append(row)

    return pd.DataFrame(rows)


def save_pyvis_html(net, html_out):
    html_out = Path(html_out)
    html_out.parent.mkdir(parents=True, exist_ok=True)
    net.save_graph(str(html_out))


def save_graphml(G, graphml_out):
    graphml_out = Path(graphml_out)
    graphml_out.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(G, str(graphml_out))


# ============================================================
# Left-side controls
# ============================================================
def network_style_controls():
    """
    Controls shown on the left side of the page.
    """
    st.subheader("Network controls")
    dg_col = DEFAULT_DG_COL
    unc_col = DEFAULT_UNC_COL
    # with st.expander("Columns", expanded=True):
    #     dg_col = st.text_input(
    #         "DG column",
    #         value=DEFAULT_DG_COL,
    #         key="rbfe_dg_col",
    #     )

        # unc_col = st.text_input(
        #     "Uncertainty column",
        #     value=DEFAULT_UNC_COL,
        #     key="rbfe_unc_col",
        # )

    with st.expander("Node style", expanded=True):
        cm1, cm2, cm3 = st.columns(3)
        with cm1:
            node_margin = st.slider(
                "Node size",
                min_value=5,
                max_value=100,
                value=28,
                step=1,
                key="rbfe_node_margin",
            )
        with cm2:
            node_font_size = st.slider(
                "Node text size",
                min_value=8,
                max_value=50,
                value=22,
                step=1,
                key="rbfe_node_font_size",
            )
        with cm3:
            node_border_width = st.slider(
                "Node border width",
                min_value=1,
                max_value=10,
                value=2,
                step=1,
                key="rbfe_node_border_width",
            )

    with st.expander("Edge style", expanded=True):
        cn11, cn21 = st.columns(2)
        with cn11:
            edge_width = st.slider(
                "Edge width",
                min_value=1,
                max_value=15,
                value=4,
                step=1,
                key="rbfe_edge_width",
            )
        with cn21:
            edge_font_size = st.slider(
                "Edge text size",
                min_value=8,
                max_value=80,
                value=20,
                step=1,
                key="rbfe_edge_font_size",
            )
        cn12, cn22 = st.columns(2)
        with cn12:
            arrow_size = st.slider(
                "Arrow size",
                min_value=0.2,
                max_value=5.0,
                value=1.5,
                step=0.1,
                key="rbfe_arrow_size",
            )
        with cn22:
            curvature = st.slider(
                "Edge curvature",
                min_value=0.0,
                max_value=1.0,
                value=0.25,
                step=0.05,
                key="rbfe_curvature",
            )

    with st.expander("Color style", expanded=True):
        cd1, cd2, cd3, cd4 = st.columns(4)
        with cd1:
            node_color = st.color_picker(
                "Node",
                "#ADD8E6",
                key="rbfe_node_color",
            )
        with cd2:
            node_border_color = st.color_picker(
                "Node border",
                "#000000",
                key="rbfe_node_border_color",
            )
        with cd3:
            positive_edge_color = st.color_picker(
                "Positive ddG",
                "#FF0000",
                key="rbfe_positive_edge_color",
            )
        with cd4:
            negative_edge_color = st.color_picker(
                "Negative ddG",
                "#0000FF",
                key="rbfe_negative_edge_color",
            )
    with st.expander("Layout", expanded=True):
        initial_layout = st.selectbox(
            "Initial layout",
            ["spring", "circular", "shell", "kamada_kawai"],
            index=0,
            key="rbfe_initial_layout",
        )
        ck1, ck2, ck3 = st.columns(3)
        with ck1:
            height_px = st.slider(
                "Canvas height",
                min_value=500,
                max_value=1600,
                value=950,
                step=50,
                key="rbfe_height_px",
            )
        with ck2:
            layout_scale = st.slider(
                "Layout spread",
                min_value=300,
                max_value=2200,
                value=1000,
                step=50,
                key="rbfe_layout_scale",
            )
        with ck3:
            rotation_deg = st.slider(
                "Rotation",
                min_value=-180,
                max_value=180,
                value=0,
                step=5,
                key="rbfe_rotation_deg",
            )

    with st.expander("Options", expanded=True):
        show_tables = st.checkbox(
            "Show summary tables",
            value=False,
            key="rbfe_show_tables",
        )

        # show_edge_titles = st.checkbox(
        #     "Show edge details on hover",
        #     value=True,
        #     key="rbfe_show_edge_titles",
        # )

        show_pyvis_controls = st.checkbox(
            "Show PyVis configure panel",
            value=False,
            key="rbfe_show_pyvis_controls",
        )

        html_filename = st.text_input(
            "HTML filename",
            value="fep_network_interactive.html",
            key="rbfe_html_filename",
        )

        # graphml_filename = st.text_input(
        #     "GraphML filename",
        #     value="fep_perturbation_network.graphml",
        #     key="rbfe_graphml_filename",
        # )

    return {
        "dg_col": dg_col,
        "unc_col": unc_col,
        "node_margin": node_margin,
        "node_font_size": node_font_size,
        "node_border_width": node_border_width,
        "node_color": node_color,
        "node_border_color": node_border_color,
        "edge_width": edge_width,
        "edge_font_size": edge_font_size,
        "arrow_size": arrow_size,
        "curvature": curvature,
        "positive_edge_color": positive_edge_color,
        "negative_edge_color": negative_edge_color,
        "initial_layout": initial_layout,
        "height_px": height_px,
        "layout_scale": layout_scale,
        "rotation_deg": rotation_deg,
        "show_tables": show_tables,
        # "show_edge_titles": show_edge_titles,
        "show_pyvis_controls": show_pyvis_controls,
        "html_filename": html_filename,
        # "graphml_filename": graphml_filename,
    }


# ============================================================
# Build PyVis network
# ============================================================
def build_pyvis_network(G, positions, controls):
    """
    Convert NetworkX graph into a draggable PyVis network.
    """
    net = Network(
        height=f"{controls['height_px']}px",
        width="100%",
        directed=True,
        notebook=False,
        bgcolor="white",
        font_color="black",
    )

    # --------------------------------------------------------
    # Nodes
    # --------------------------------------------------------
    for node in G.nodes():
        x = positions[node]["x"]
        y = positions[node]["y"]

        net.add_node(
            node,
            label=node,
            title=node,
            x=x,
            y=y,
            fixed=False,
            physics=False,
            shape="circle",
            margin=int(controls["node_margin"]),
            shapeProperties={
                "borderRadius": 100,
            },
            color={
                "background": controls["node_color"],
                "border": controls["node_border_color"],
                "highlight": {
                    "background": controls["node_color"],
                    "border": controls["node_border_color"],
                },
                "hover": {
                    "background": controls["node_color"],
                    "border": controls["node_border_color"],
                },
            },
            borderWidth=int(controls["node_border_width"]),
            font={
                "size": int(controls["node_font_size"]),
                "face": "Arial",
                "color": "black",
                "align": "center",
                "vadjust": 0,
            },
        )

    # --------------------------------------------------------
    # Edges
    # --------------------------------------------------------
    for u, v, d in G.edges(data=True):
        if "ddG" not in d:
            continue

        ddg = float(d["ddG"])
        ddg_unc = float(d["ddG_unc"])

        color = (
            controls["positive_edge_color"]
            if ddg > 0
            else controls["negative_edge_color"]
        )

        label = f"{ddg:.2f} ± {ddg_unc:.2f}"

        # if controls["show_edge_titles"]:
        #     title = (
        #         f"{u} → {v}<br>"
        #         f"ddG = {ddg:.3f} ± {ddg_unc:.3f} kcal/mol"
        #     )

        #     if "complex_DG" in d and "complex_unc" in d:
        #         title += (
        #             f"<br>complex DG = "
        #             f"{d['complex_DG']:.3f} ± {d['complex_unc']:.3f}"
        #         )

        #     if "solvent_DG" in d and "solvent_unc" in d:
        #         title += (
        #             f"<br>solvent DG = "
        #             f"{d['solvent_DG']:.3f} ± {d['solvent_unc']:.3f}"
        #         )
        # else:
        #     title = ""

        net.add_edge(
            u,
            v,
            label=label,
            # title=title,
            color=color,
            width=int(controls["edge_width"]),
            arrows={
                "to": {
                    "enabled": True,
                    "scaleFactor": float(controls["arrow_size"]),
                    "type": "arrow",
                }
            },
            font={
                # "size": int(controls["edge_font_size"]),
                "face": "Arial",
                "color": color,
                "background": "white",
                "strokeWidth": 0,
                "align": "middle",
            },
            smooth={
                "enabled": True,
                "type": "curvedCW",
                "roundness": float(controls["curvature"]),
            },
        )

    # --------------------------------------------------------
    # Global options
    # --------------------------------------------------------
    options = {
        "autoResize": True,
        "layout": {
            "improvedLayout": False,
        },
        "interaction": {
            "dragNodes": True,
            "dragView": True,
            "zoomView": True,
            "hover": True,
            "navigationButtons": True,
            "keyboard": True,
        },
        "physics": {
            "enabled": False,
        },
        "nodes": {
            "shape": "box",
            "margin": int(controls["node_margin"]),
            "font": {
                "size": int(controls["node_font_size"]),
                "face": "Arial",
                "color": "black",
                "align": "center",
                "vadjust": 0,
                "multi": False,
            },
        },
        "edges": {
            "font": {
                "size": int(controls["edge_font_size"]),
                "face": "Arial",
                "align": "middle",
                "strokeWidth": 0,
            },
            "smooth": {
                "enabled": True,
                "type": "curvedCW",
                "roundness": float(controls["curvature"]),
            },
        },
    }

    if controls["show_pyvis_controls"]:
        options["configure"] = {
            "enabled": True,
            "filter": [
                "nodes",
                "edges",
                "layout",
                "physics",
            ],
            "showButton": True,
        }

    net.set_options(json.dumps(options))

    return net


# ============================================================
# Main tab function
# ============================================================
def tab2_design(analysis_tsv_picker, workdir):
    st.header("RBFE Network")

    # --------------------------------------------------------
    # Select TSV file
    # --------------------------------------------------------
    selected_tsv_file = analysis_tsv_picker(
        label="Input analysis TSV file",
        start_dir=str(workdir),
        allowed_extensions=(".tsv",),
        key_prefix="raw_tsv",
    )

    st.divider()

    if selected_tsv_file is None:
        st.info(
            "Click **List**, browse to your analysis TSV file, "
            "then click **Use TSV file**."
        )
        return

    tsv_file = Path(selected_tsv_file)

    if not tsv_file.exists():
        st.error(f"Cannot find TSV file:\n\n{tsv_file}")
        return

    try:
        df = pd.read_csv(str(tsv_file), sep="\t")
    except Exception as e:
        st.error(f"Could not read TSV file:\n\n{e}")
        return

    st.success(f"Loaded TSV file: {tsv_file}")

    with st.expander("TSV preview", expanded=False):
        st.dataframe(df.head(), use_container_width=True)

    # --------------------------------------------------------
    # Left controls, right network
    # --------------------------------------------------------
    left_col, right_col = st.columns([1, 3], gap="large")

    with left_col:
        controls = network_style_controls()

    # --------------------------------------------------------
    # Validate input columns
    # --------------------------------------------------------
    required_cols = [
        "leg",
        "ligand_i",
        "ligand_j",
        controls["dg_col"],
        controls["unc_col"],
    ]

    missing_cols = [c for c in required_cols if c not in df.columns]

    if missing_cols:
        with right_col:
            st.error(f"Missing required columns:\n\n{missing_cols}")
            st.write("Available columns:")
            st.write(list(df.columns))
        return

    # --------------------------------------------------------
    # Build graph
    # --------------------------------------------------------
    try:
        G, summary = build_fep_graph(
            df,
            controls["dg_col"],
            controls["unc_col"],
        )
    except Exception as e:
        with right_col:
            st.error(f"Could not build FEP graph:\n\n{e}")
        return

    edge_df = make_edge_table(G)

    positions = make_initial_positions(
        G,
        controls["initial_layout"],
        controls["layout_scale"],
        rotation_deg=controls["rotation_deg"],
    )

    net = build_pyvis_network(
        G=G,
        positions=positions,
        controls=controls,
    )

    # --------------------------------------------------------
    # Output paths
    # --------------------------------------------------------
    if tsv_file.parent.name == "analysis":
        analysis_dir = tsv_file.parent
    else:
        analysis_dir = tsv_file.parent / "analysis"

    html_out = analysis_dir / controls["html_filename"]
    # graphml_out = analysis_dir / controls["graphml_filename"]

    # --------------------------------------------------------
    # Right side: buttons and graph
    # --------------------------------------------------------
    with right_col:
        st.subheader("Interactive RBFE network")

        b1, b2, b3 = st.columns(3)

        with b1:
            if st.button("Save interactive HTML", key="rbfe_save_html"):
                try:
                    save_pyvis_html(net, html_out)
                    st.success(f"Saved HTML to:\n\n{html_out}")
                except Exception as e:
                    st.error(f"Could not save HTML:\n\n{e}")

        with b2:
            if st.button("Save GraphML", key="rbfe_save_graphml"):
                try:
                    save_graphml(G, graphml_out)
                    st.success(f"Saved GraphML to:\n\n{graphml_out}")
                except Exception as e:
                    st.error(f"Could not save GraphML:\n\n{e}")

        with b3:
            st.info(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

        with tempfile.NamedTemporaryFile(
            delete=False,
            suffix=f"_{time.time_ns()}.html",
        ) as tmp:
            net.save_graph(tmp.name)
            html = Path(tmp.name).read_text(encoding="utf-8")

        position_storage_key = f"rbfe_positions_{tsv_file.stem}"

        html = add_persist_node_positions_js(
            html,
            storage_key=position_storage_key,
        )

        html = add_download_png_button(html)

        components.html(
            html,
            height=int(controls["height_px"]) + 120,
            scrolling=True,
        )

    # --------------------------------------------------------
    # Optional tables
    # --------------------------------------------------------
    if controls["show_tables"]:
        st.divider()

        st.subheader("Per-leg averaged edges")
        st.dataframe(summary, use_container_width=True)

        st.subheader("Computed ddG edges")
        st.dataframe(edge_df, use_container_width=True)
