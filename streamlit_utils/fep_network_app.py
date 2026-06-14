import os
import json
import tempfile

import numpy as np
import pandas as pd
import networkx as nx
import streamlit as st
import streamlit.components.v1 as components
from pyvis.network import Network


tsv_file = st.text_input(
    "Input TSV file",
    "/data/liuy48/openfe_m/simulation/rbfe/PDE4B/A33/raw_new.tsv",
)


html_out = st.text_input(
    "Output HTML file",
    "/data/liuy48/openfe_m/simulation/rbfe/PDE4B/A33/analysis/fep_network_interactive.html",
)

st.header("Column names")

dg_col = st.text_input(
    "DG column",
    "DG(i->j) (kcal/mol)",
)

unc_col = st.text_input(
    "Uncertainty column",
    "MBAR uncertainty (kcal/mol)",
)


# ============================================================
# Visual controls
# ============================================================
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.header("Node style")
    node_margin = st.slider(
        "Node padding / size",
        min_value=5,
        max_value=100,
        value=28,
        step=1,
    )

    node_font_size = st.slider(
        "Node text size",
        min_value=8,
        max_value=50,
        value=22,
        step=1,
    )

    node_border_width = st.slider(
        "Node border width",
        min_value=1,
        max_value=10,
        value=2,
        step=1,
    )

with c2:
    st.header("Edge style")

    edge_width = st.slider(
        "Edge width",
        min_value=1,
        max_value=15,
        value=4,
        step=1,
    )

    edge_font_size = st.slider(
        "Edge label text size",
        min_value=8,
        max_value=80,
        value=20,
        step=1,
    )

    arrow_size = st.slider(
        "Arrow size",
        min_value=0.2,
        max_value=5.0,
        value=1.5,
        step=0.1,
    )

    curvature = st.slider(
        "Edge curvature",
        min_value=0.0,
        max_value=1.0,
        value=0.25,
        step=0.05,
    )
with c3:
    st.header("Canvas / layout")

    initial_layout = st.selectbox(
        "Initial layout",
        ["spring", "circular", "shell", "kamada_kawai"],
        index=0,
    )
    height_px = st.slider(
        "Canvas height",
        min_value=500,
        max_value=1600,
        value=950,
        step=50,
    )

    layout_scale = st.slider(
        "Layout spread",
        min_value=300,
        max_value=2200,
        value=1000,
        step=50,
    )

    rotation_deg = st.slider(
        "Rotate initial layout degrees",
        min_value=-180,
        max_value=180,
        value=0,
        step=5,
    )
with c4:
    st.header("Options")

    show_buttons = st.checkbox(
        "Show PyVis controls",
        value=True,
    )

    show_tables = st.checkbox(
        "Show summary tables",
        value=True,
    )

    show_edge_titles = st.checkbox(
        "Show edge details on hover",
        value=True,
    )


# ============================================================
# Helper functions
# ============================================================
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

    html = html.replace("</body>", download_js + "\n</body>")
    return html

def make_initial_positions(G, layout_name, scale, rotation_deg=0):
    """
    Generate initial node positions for PyVis.

    The rotation only affects the starting layout.
    After rendering, you can still drag nodes manually.
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

        # Rotate around origin.
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


def build_pyvis_network(
    G,
    positions,
    height_px,
    node_margin,
    node_font_size,
    node_border_width,
    edge_width,
    edge_font_size,
    arrow_size,
    curvature,
    show_buttons,
    show_edge_titles,
):
    """
    Convert NetworkX graph into a draggable PyVis network.
    """
    net = Network(
        height=f"{height_px}px",
        width="100%",
        directed=True,
        notebook=False,
        bgcolor="white",
        font_color="black",
    )
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        node_color = st.color_picker("Node color", "#ADD8E6")
    with c2:
        node_border_color = st.color_picker("Node border color", "#000000")
    with c3: 
        positive_edge_color = st.color_picker("Positive ddG edge color", "#FF0000")
    with c4: 
        negative_edge_color = st.color_picker("Negative ddG edge color", "#0000FF")
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

            # The node can be dragged manually.
            fixed=False,

            # Important: disables physics per node.
            # This prevents other nodes from moving when one node is dragged.
            physics=False,

            # Circle with internal label.
            shape="circle",

            # For circle nodes with text inside, margin controls visible node size.
            margin=int(node_margin),

            color={
                "background": node_color,
                "border": node_border_color,
            },

            borderWidth=int(node_border_width),

            font={
                "size": int(node_font_size),
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

        
        color = positive_edge_color if ddg > 0 else negative_edge_color
        label = f"{ddg:.2f} ± {ddg_unc:.2f}"

        if show_edge_titles:
            title = (
                f"{u} → {v}<br>"
                f"ddG = {ddg:.3f} ± {ddg_unc:.3f} kcal/mol"
            )

            if "complex_DG" in d and "complex_unc" in d:
                title += (
                    f"<br>complex DG = "
                    f"{d['complex_DG']:.3f} ± {d['complex_unc']:.3f}"
                )

            if "solvent_DG" in d and "solvent_unc" in d:
                title += (
                    f"<br>solvent DG = "
                    f"{d['solvent_DG']:.3f} ± {d['solvent_unc']:.3f}"
                )
        else:
            title = ""

        net.add_edge(
            u,
            v,
            label=label,
            title=title,
            color=color,
            width=int(edge_width),

            # Arrow size is controlled here.
            arrows={
                "to": {
                    "enabled": True,
                    "scaleFactor": float(arrow_size),
                    "type": "arrow",
                }
            },

            font={
                "size": int(edge_font_size),
                "face": "Arial",
                "color": color,
                "background": "white",
                "strokeWidth": 0,
                "align": "middle",
            },

            smooth={
                "enabled": True,
                "type": "curvedCW",
                "roundness": curvature,
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
            "margin": int(node_margin),
            "font": {
                "size": int(node_font_size),
                "face": "Arial",
                "color": "black",
                "align": "center",
                "vadjust": 0,
                "multi": False,
            },
        },
        "edges": {
            "font": {
                "size": int(edge_font_size),
                "face": "Arial",
                "align": "middle",
                "strokeWidth": 0,
            },
            "smooth": {
                "enabled": True,
                "type": "curvedCW",
                "roundness": float(curvature),
            },
        },
    }

    if show_buttons:
        options["configure"] = {
            "enabled": True,
            "filter": [
                "nodes",
                "edges",
                "layout",
                "physics",
            ],
        }

    net.set_options(json.dumps(options))

    return net


# def save_graphml(G, out_graphml):
#     out_dir = os.path.dirname(out_graphml)
#     if out_dir:
#         os.makedirs(out_dir, exist_ok=True)

#     nx.write_graphml(G, out_graphml)


def save_pyvis_html(net, html_out):
    out_dir = os.path.dirname(html_out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    net.save_graph(html_out)


# ============================================================
# Main app
# ============================================================
if not os.path.exists(tsv_file):
    st.error(f"Cannot find TSV file:\n\n{tsv_file}")
    st.stop()

try:
    df = pd.read_csv(tsv_file, sep="\t")
except Exception as e:
    st.error(f"Could not read TSV file:\n\n{e}")
    st.stop()

required_cols = [
    "leg",
    "ligand_i",
    "ligand_j",
    dg_col,
    unc_col,
]

missing_cols = [c for c in required_cols if c not in df.columns]

if missing_cols:
    st.error(f"Missing required columns:\n\n{missing_cols}")
    st.write("Available columns:")
    st.write(list(df.columns))
    st.stop()

try:
    G, summary = build_fep_graph(df, dg_col, unc_col)
except Exception as e:
    st.error(f"Could not build FEP graph:\n\n{e}")
    st.stop()

edge_df = make_edge_table(G)

positions = make_initial_positions(
    G,
    initial_layout,
    layout_scale,
    rotation_deg=rotation_deg,
)

# try:
#     save_graphml(G, out_graphml)
# except Exception as e:
#     st.warning(f"Could not save GraphML automatically:\n\n{e}")

net = build_pyvis_network(
    G=G,
    positions=positions,
    height_px=height_px,
    node_margin=node_margin,
    node_font_size=node_font_size,
    node_border_width=node_border_width,
    edge_width=edge_width,
    edge_font_size=edge_font_size,
    arrow_size=arrow_size,
    curvature=curvature,
    show_buttons=show_buttons,
    show_edge_titles=show_edge_titles,
)


# ============================================================
# Buttons
# ============================================================
col1, col2, col3 = st.columns(3)

with col1:
    if st.button("Save interactive HTML"):
        try:
            save_pyvis_html(net, html_out)
            st.success(f"Saved HTML to:\n\n{html_out}")
        except Exception as e:
            st.error(f"Could not save HTML:\n\n{e}")


with col3:
    st.info(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")


# ============================================================
# Render interactive graph
# ============================================================
with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp:
    net.save_graph(tmp.name)
    html = open(tmp.name, "r", encoding="utf-8").read()

html = add_download_png_button(html)

components.html(
    html,
    height=height_px + 120,
    scrolling=True,
)


# ============================================================
# Tables
# ============================================================
if show_tables:
    st.subheader("Per-leg averaged edges")
    st.dataframe(summary, use_container_width=True)

    st.subheader("Computed ddG edges")
    st.dataframe(edge_df, use_container_width=True)
