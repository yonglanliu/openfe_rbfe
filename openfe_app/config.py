from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @property
    def inputs(self) -> Path:
        return self.root / "inputs"

    @property
    def planned(self) -> Path:
        return self.root / "json"   # output from plan-rbfe-network

    @property
    def results(self) -> Path:
        return self.root / "results"        # output jsons from quickrun

    @property
    def analysis(self) -> Path:
        return self.root / "analysis"       # gathered tables/plots

    def ensure(self) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        self.inputs.mkdir(exist_ok=True)
        self.planned.mkdir(exist_ok=True)
        self.results.mkdir(exist_ok=True)
        self.analysis.mkdir(exist_ok=True)
