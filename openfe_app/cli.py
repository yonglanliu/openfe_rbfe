from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import subprocess
import shlex
import os


@dataclass
class CmdResult:
    returncode: int
    stdout: str
    stderr: str
    cmd: List[str]


def run_cmd(cmd: List[str], cwd: Optional[Path] = None, env: Optional[dict] = None) -> CmdResult:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return CmdResult(proc.returncode, proc.stdout, proc.stderr, cmd)


def openfe_plan_rbfe_network(
    ligands_sdf: Path,
    protein_pdb: Path,
    out_dir: Path,
    settings_yaml: Optional[Path] = None,
) -> CmdResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "openfe", "plan-rbfe-network",
        "-M", str(ligands_sdf),
        "-p", str(protein_pdb),
        "-o", str(out_dir),
    ]
    if settings_yaml:
        cmd.extend(["-s", str(settings_yaml)])
    return run_cmd(cmd)


def list_transformations(planned_dir: Path) -> List[Path]:
    # OpenFE CLI writes transformation json files into output dir (as shown in docs)
    return sorted(planned_dir.glob("*.json"))


def classify_leg(json_path: Path) -> str:
    name = json_path.name.lower()
    if "complex" in name:
        return "complex"
    if "solvent" in name:
        return "solvent"
    if "vacuum" in name:
        return "vacuum"
    return "unknown"


def openfe_quickrun(
    transformation_json: Path,
    out_result_json: Path,
    workdir: Path,
) -> CmdResult:
    workdir.mkdir(parents=True, exist_ok=True)
    out_result_json.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "openfe", "quickrun",
        str(transformation_json),
        "-o", str(out_result_json),
        "-d", str(workdir),
    ]
    return run_cmd(cmd)


def openfe_gather(results_dir: Path, out_tsv: Path, report: str = "dg") -> CmdResult:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    # Docs show running gather from within working directory, but also accepts path argument
    cmd = ["openfe", "gather", str(results_dir), "--report", report, "-o", str(out_tsv)]
    return run_cmd(cmd)
