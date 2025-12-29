"""
openfe_app/view3d.py

Utilities + Streamlit/py3Dmol viewer for:
- Extracting ligands from PDB (HETATM grouped by residue)
- Loading ligands from SDF as RDKit mols
- Aligning ligands to a reference (MCS-based)
- Rendering ONE interactive 3D window with configurable protein/ligand representations

Dependencies:
    pip install streamlit py3Dmol
    conda install -c conda-forge rdkit
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Tuple, Dict
import pandas as pd
import streamlit as st
import py3Dmol
from streamlit.components.v1 import html


# ============================
#  A) Chemistry / IO utilities
# ============================

@dataclass
class PDBLigandBlock:
    resname: str
    chain: str
    resid: str
    pdb_block: str


def _is_water_resname(resn: str) -> bool:
    return resn.upper() in {"HOH", "WAT", "H2O", "SPC"}


def _is_common_ion(resn: str) -> bool:
    return resn.upper() in {
        "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "BR", "I",
        "CS", "SR", "CO", "NI", "CD", "HG"
    }


def extract_ligands_from_pdb(pdb_path: Union[str, Path]) -> List[PDBLigandBlock]:
    """
    Extract non-protein residues from a PDB file as per-residue PDB blocks (HETATM groups),
    excluding waters and common ions.

    Returns list sorted by decreasing atom count (largest first).
    """
    pdb_path = Path(pdb_path)
    lines = pdb_path.read_text(errors="ignore").splitlines()

    groups: Dict[Tuple[str, str, str], List[str]] = {}
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        
        resname = line[17:20].strip()
        chain = (line[21].strip() or " ")
        resseq = line[22:26].strip()
        icode = line[26].strip()
        resid = f"{resseq}{icode}".strip() if icode else resseq

        if _is_water_resname(resname) or _is_common_ion(resname):
            continue

        key = (resname, chain, resid)
        groups.setdefault(key, []).append(line)
    #st.info(groups)

    ligs: List[PDBLigandBlock] = []
    for (resname, chain, resid), het_lines in groups.items():
        if not (resname == "NME" or resname =="ACE"):
            block = "\n".join(het_lines) + "\nEND\n"
            ligs.append(PDBLigandBlock(resname=resname, chain=chain, resid=resid, pdb_block=block))

    ligs.sort(key=lambda x: x.pdb_block.count("\n"), reverse=True)
    return ligs


def rdkit_mol_from_pdb_block(pdb_block: str):
    """
    Convert a small-molecule PDB block to an RDKit Mol.
    Bond orders in PDB are ambiguous; RDKit may fail for some ligands.
    """
    try:
        from rdkit import Chem
    except Exception as e:
        raise RuntimeError("RDKit required. Install: conda install -c conda-forge rdkit") from e

    m = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        try:
            Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
    return m


def rdkit_mols_from_sdf(sdf_path: Union[str, Path]) -> List["Chem.Mol"]:
    """Read molecules from an SDF file path into a list of RDKit Mol."""
    try:
        from rdkit import Chem
    except Exception as e:
        raise RuntimeError("RDKit required. Install: conda install -c conda-forge rdkit") from e

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    return [m for m in suppl if m is not None]

# Extract name of ligands
def mol_name(m, fallback: str = "lig") -> str:
    """Get a display name for an RDKit Mol."""
    try:
        if m.HasProp("_Name"):
            nm = m.GetProp("_Name").strip()
            if nm:
                return nm
    except Exception:
        pass
    return fallback


def _ensure_3d_conformer(m):
    """Ensure mol has a 3D conformer; if not, embed one."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m2 = Chem.Mol(m)
    if m2.GetNumConformers() == 0:
        m2 = Chem.AddHs(m2, addCoords=True)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xC0FFEE
        ok = AllChem.EmbedMolecule(m2, ps)
        if ok != 0:
            ok = AllChem.EmbedMolecule(m2, AllChem.ETKDG())
        if ok == 0:
            try:
                AllChem.UFFOptimizeMolecule(m2, maxIters=200)
            except Exception:
                pass
    return m2


def align_to_reference(mols: List["Chem.Mol"], ref: "Chem.Mol") -> List["Chem.Mol"]:
    """
    Align mols to ref using MCS-based atom mapping.
    Returns aligned copies in same order.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS

    if not mols:
        return []
    ref3d = _ensure_3d_conformer(ref)

    aligned: List[Chem.Mol] = []
    for m in mols:
        m3d = _ensure_3d_conformer(m)

        if ref3d.GetNumConformers() == 0 or m3d.GetNumConformers() == 0:
            aligned.append(m3d)
            continue

        try:
            mcs = rdFMCS.FindMCS(
                [ref3d, m3d],
                timeout=20,
                ringMatchesRingOnly=True,
                completeRingsOnly=True,
                matchValences=True,
            )
            if not mcs.smartsString:
                aligned.append(m3d)
                continue

            patt = Chem.MolFromSmarts(mcs.smartsString)
            ref_match = ref3d.GetSubstructMatch(patt)
            mol_match = m3d.GetSubstructMatch(patt)

            if not ref_match or not mol_match or len(ref_match) != len(mol_match):
                aligned.append(m3d)
                continue

            atom_map = list(zip(mol_match, ref_match))
            AllChem.AlignMol(m3d, ref3d, atomMap=atom_map)
            aligned.append(m3d)
        except Exception:
            aligned.append(m3d)

    return aligned

# ============================
#  B) Find pocket residues
# ============================
import math
from collections import defaultdict

def _pdb_atoms_by_residue(pdb_text: str):
    # returns dict[(chain, resi)] -> list of (x,y,z)
    res_atoms = defaultdict(list)
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        # skip hetero for protein pocket calc (waters/ions/ligands in pdb)
        if line.startswith("HETATM"):
            continue

        chain = line[21].strip() or " "
        resi = line[22:26].strip()
        if not resi:
            continue
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        res_atoms[(chain, resi)].append((x, y, z))
    return res_atoms

def _sdf_atom_coords(sdf_text: str):
    # very small SDF parser: reads counts line then N atom lines
    lines = sdf_text.splitlines()
    if len(lines) < 4:
        return []
    counts = lines[3]
    try:
        n_atoms = int(counts[0:3])
    except:
        return []
    coords = []
    start = 4
    for i in range(n_atoms):
        parts = lines[start + i].split()
        if len(parts) < 3:
            continue
        coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return coords

def _find_pocket_residues(pdb_text: str, sdf_texts: list[str], radius: float):
    r2 = radius * radius
    res_atoms = _pdb_atoms_by_residue(pdb_text)

    lig_coords = []
    for sdf in sdf_texts:
        lig_coords.extend(_sdf_atom_coords(sdf))

    pocket = set()
    if not lig_coords:
        return pocket

    # brute force (fine for typical sizes)
    for (chain, resi), atoms in res_atoms.items():
        for ax, ay, az in atoms:
            for lx, ly, lz in lig_coords:
                dx = ax - lx; dy = ay - ly; dz = az - lz
                if dx*dx + dy*dy + dz*dz <= r2:
                    pocket.add((chain, resi))
                    break
            if (chain, resi) in pocket:
                break
    return pocket


# ============================
#  C) ONE-window 3D viewer
# ============================

def _show_py3dmol(view: py3Dmol.view, height: int = 680):
    html(view._make_html(), height=height, width=None)


def _as_sdf_texts(ligands: List[Union[str, "Chem.Mol"]]) -> List[str]:
    if not ligands:
        return []
    if isinstance(ligands[0], str):
        return [x for x in ligands if isinstance(x, str)]
    try:
        from rdkit import Chem  # type: ignore
    except Exception as e:
        raise RuntimeError("RDKit required for Mol inputs.") from e
    return [Chem.MolToMolBlock(m) for m in ligands if m is not None]


def render_protein_ligands_single_view(
    protein_pdb_text: str,
    aligned_ligands: List[Union[str, "Chem.Mol"]],
    ligand_names: Optional[List[str]] = None,
    key_prefix: str = "single3d",
    width: int = 1100,
    height: int = 680,
):
    """
    ONE viewer window:
      - protein representation options: cartoon/lines/sticks/sphere/surface/cartoon+surface
      - ligand representation options: sticks/lines/sphere/surface
      - show ALL aligned ligands by default (multiselect if you want to hide some)
      - pocket highlight (protein residues within radius of any shown ligand)
    """
    if not protein_pdb_text or not protein_pdb_text.strip():
        st.info("Upload protein PDB first.")
        return

    ligands_sdf = _as_sdf_texts(aligned_ligands)
    n_ligs = len(ligands_sdf)
    if ligand_names is None or (n_ligs and len(ligand_names) != n_ligs):
        ligand_names = [f"ligand_{i+1}" for i in range(n_ligs)]

    # ---- Controls (protein + ligand) ----
    with st.expander("3D display options", expanded=True):
        c1, c2, c3 = st.columns([1, 1, 1])

        with c1:
            protein_repr = st.selectbox(
                "Protein representation",
                ["cartoon", "lines", "sticks", "sphere", "surface", "cartoon+surface"],
                index=0,
                key=f"{key_prefix}_protein_repr",
            )
            protein_color = st.selectbox(
                "Protein color",
                ["spectrum", "white", "gray", "lightgray", "slateblue"],
                index=0,
                key=f"{key_prefix}_protein_color",
            )
            surface_opacity = 0.65
            if "surface" in protein_repr:
                surface_opacity = st.slider(
                    "Protein surface opacity",
                    0.05, 1.0, 0.65, 0.05,
                    key=f"{key_prefix}_surf_opacity",
                )

        with c2:
            ligand_repr = st.selectbox(
                "Ligand representation",
                ["sticks", "lines", "sphere", "surface"],
                index=0,
                key=f"{key_prefix}_ligand_repr",
            )
            ligand_size = st.slider(
                "Ligand size",
                0.10, 1.00, 0.25, 0.05,
                key=f"{key_prefix}_ligand_size",
            )
            ligand_color_mode = st.selectbox(
                "Ligand coloring",
                ["auto (different each)", "single color"],
                index=0,
                key=f"{key_prefix}_lig_color_mode",
            )
            single_lig_color = "green"
            if ligand_color_mode == "single color":
                single_lig_color = st.selectbox(
                    "Single ligand color",
                    ["green", "cyan", "magenta", "orange", "yellow", "red", "blue", "white", "gray"],
                    index=0,
                    key=f"{key_prefix}_lig_single_color",
                )

        with c3:
            show_pocket = st.toggle(
                "Highlight pocket residues near ligands",
                value=False,
                key=f"{key_prefix}_show_pocket",
            )
            pocket_radius = 5.0
            if show_pocket:
                pocket_radius = st.slider(
                    "Pocket radius (Å)",
                    2.0, 12.0, 5.0, 0.5,
                    key=f"{key_prefix}_pocket_radius",
                )

            # If many ligands, allow hiding some; default shows ALL (what you asked)
            if n_ligs <= 120:
                default_sel = ligand_names
            else:
                default_sel = ligand_names[:120]

            selected_ligands = st.multiselect(
                "Ligands to display (default = all)",
                options=ligand_names,
                default=default_sel,
                key=f"{key_prefix}_selected_ligs",
            )
            

    # map selection -> indices
    idx_by_name = {nm: i for i, nm in enumerate(ligand_names)}
    selected_indices = [idx_by_name[nm] for nm in selected_ligands if nm in idx_by_name]

    # ---- Build viewer ----
    view = py3Dmol.view(width=width, height=height)
    view.addModel(protein_pdb_text, "pdb")
    view.setStyle({}, {})  # reset

   # --- Protein styling (robust): always model 0 ---
    protein_sel = {"model": 0}  # model 0 is the protein PDB we added first

    # Optional: hide hetero atoms from the protein model (waters/ions/ligand in PDB)
    # This prevents the PDB ligand from being shown twice (once from PDB, once from aligned SDF)
    view.setStyle({"model": 0, "hetflag": True}, {})  # clear het styling
    view.setStyle({"model": 0, "hetflag": False}, {})  # clear protein styling

    if protein_color == "spectrum":
        pcolor = {"color": "spectrum"}
    else:
        pcolor = {"color": protein_color}

    if protein_repr == "cartoon":
        view.setStyle({"model": 0, "hetflag": False}, {"cartoon": pcolor})

    elif protein_repr == "lines":
        view.setStyle({"model": 0, "hetflag": False}, {"line": pcolor})

    elif protein_repr == "sticks":
        view.setStyle({"model": 0, "hetflag": False}, {"stick": {"colorscheme": "elem"}})

    elif protein_repr == "sphere":
        view.setStyle({"model": 0, "hetflag": False}, {"sphere": {"scale": 0.3}})

    elif protein_repr == "surface":
        # show cartoon under surface for orientation
        view.addSurface(py3Dmol.VDW, {"opacity": surface_opacity, **pcolor}, {"model": 0, "hetflag": False})

    elif protein_repr == "cartoon+surface":
        view.setStyle({"model": 0, "hetflag": False}, {"cartoon": pcolor})
        view.addSurface(py3Dmol.VDW, {"opacity": surface_opacity, **pcolor}, {"model": 0, "hetflag": False})


    # Ligands: add as separate models so they can have independent colors
    palette = ["green", "cyan", "magenta", "orange", "yellow", "red", "blue", "lime", "pink", "purple"]
    model_id = 1  # protein is model 0
    ligand_model_ids: List[int] = []

    for j, lig_i in enumerate(selected_indices):
        view.addModel(ligands_sdf[lig_i], "sdf")
        ligand_model_ids.append(model_id)

        lig_color = palette[j % len(palette)] if ligand_color_mode == "auto (different each)" else single_lig_color
        sel = {"model": model_id}

        if ligand_repr == "sticks":
            view.setStyle(sel, {"stick": {"color": lig_color, "radius": ligand_size}})
        elif ligand_repr == "lines":
            view.setStyle(sel, {"line": {"color": lig_color}})
        elif ligand_repr == "sphere":
            view.setStyle(sel, {"sphere": {"color": lig_color, "scale": max(0.1, ligand_size)}})
        elif ligand_repr == "surface":
            # show stick underneath + surface over
            view.setStyle(sel, {"stick": {"color": lig_color, "radius": max(0.12, ligand_size)}})
            view.addSurface(py3Dmol.VDW, {"opacity": 0.85, "color": lig_color}, sel)

        model_id += 1

    # Pocket highlight: protein residues within pocket_radius Å of ANY ligand atoms
    # Pocket highlight: protein residues within pocket_radius Å of ANY ligand atoms
    if show_pocket and selected_indices:
        # compute using the *same* ligands you're displaying
        sdf_texts = [ligands_sdf[i] for i in selected_indices]
        pocket = _find_pocket_residues(protein_pdb_text, sdf_texts, float(pocket_radius))

        # ---- Build resname lookup inline (no new function) ----
        resname_by_chain_resi = {}
        for line in protein_pdb_text.splitlines():
            if not line.startswith("ATOM  "):
                continue
            chain = (line[21].strip() or " ")
            resi = line[22:26].strip()
            resname = line[17:20].strip()
            if resi:
                resname_by_chain_resi[(chain, resi)] = resname

        # ---- Streamlit block: pocket residue list ----
        pocket_rows = []
        for chain, resi in sorted(pocket, key=lambda x: (x[0], int(x[1]))):
            pocket_rows.append({
                "Chain": chain,
                "Resid": int(resi),
                "Resname": resname_by_chain_resi.get((chain, resi), "UNK"),
            })


        # ---- Highlight pocket sticks by residue selection ----
        for chain, resi in pocket:
            view.addStyle(
                {"model": 0, "hetflag": False, "chain": chain, "resi": int(resi)},
                {"stick": {"color": "yellow", "radius": 0.25}}
            )


    view.setBackgroundColor("white")
    view.zoomTo()

    # ---- Layout: structure (left) | pocket table (right) ----
    col_view, col_table = st.columns([3, 1], gap="medium")

    with col_view:
        _show_py3dmol(view, height=height)

    with col_table:
        if show_pocket and selected_indices:
            if pocket_rows:
                st.markdown("### Pocket residues")
                st.dataframe(
                    pd.DataFrame(pocket_rows),
                    use_container_width=True,
                    height=height - 120,
                )
            else:
                st.info("No pocket residues found.")
