from __future__ import annotations
from pathlib import Path
from typing import Optional
import shutil


def save_upload(uploaded_file, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return dest


def safe_copy(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def find_first(path: Path, patterns: list[str]) -> Optional[Path]:
    for pat in patterns:
        hits = sorted(path.glob(pat))
        if hits:
            return hits[0]
    return None
