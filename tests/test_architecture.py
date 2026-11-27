"""Vérifie les règles de dépendance entre modules décrites dans docs/ARCHITECTURE.md.

Ce test est la version exécutable de la matrice de dépendances : si un module
importe un module qu'il n'a pas le droit de connaître, le test échoue et
indique le fichier fautif. Ajouter un nouveau module impose de le déclarer ici
et dans docs/ARCHITECTURE.md (section 3).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "pepforge"

# Module → modules internes qu'il a le droit d'importer.
# Ordre de lecture : des couches basses (sans dépendance) vers les couches hautes.
ALLOWED: dict[str, set[str]] = {
    "exceptions": set(),
    "core": {"exceptions"},
    "config": {"exceptions"},
    "tools": {"exceptions"},
    "residues": {"core", "exceptions"},
    "checks": {"core", "config", "exceptions"},
    "report": {"core", "config", "checks", "exceptions"},
    "io": {"core", "residues", "report", "exceptions"},
    "build": {"core", "residues", "report", "tools", "exceptions"},
    "forcefields": {"core", "config", "residues", "report", "tools", "exceptions"},
    "engines": {"core", "config", "forcefields", "report", "tools", "exceptions"},
    "pipeline": {
        "core",
        "config",
        "residues",
        "checks",
        "report",
        "io",
        "build",
        "forcefields",
        "engines",
        "tools",
        "exceptions",
    },
    "cli": {"__root__", "pipeline", "config", "io", "core", "checks", "report", "exceptions"},
}


def _unit_of(path: Path) -> str:
    relative = path.relative_to(SRC)
    if len(relative.parts) == 1:
        return "__root__" if relative.stem == "__init__" else relative.stem
    return relative.parts[0]


def _internal_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    targets: set[str] = set()
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names = [node.module]
        elif isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        for name in names:
            parts = name.split(".")
            if parts[0] != "pepforge":
                continue
            targets.add(parts[1] if len(parts) > 1 else "__root__")
    return targets


def _source_files() -> list[Path]:
    return sorted(p for p in SRC.rglob("*.py") if _unit_of(p) != "__root__")


def test_every_module_is_declared() -> None:
    units = {_unit_of(p) for p in _source_files()}
    undeclared = units - set(ALLOWED)
    assert not undeclared, f"Modules absents de la matrice de dépendances : {sorted(undeclared)}"


@pytest.mark.parametrize("path", _source_files(), ids=lambda p: str(p.relative_to(SRC)))
def test_dependencies_respect_the_matrix(path: Path) -> None:
    unit = _unit_of(path)
    forbidden = _internal_imports(path) - ALLOWED[unit] - {unit}
    assert not forbidden, f"{path.relative_to(SRC)} ({unit}) importe {sorted(forbidden)}"


def test_no_relative_imports() -> None:
    offenders = []
    for path in _source_files():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        if any(isinstance(n, ast.ImportFrom) and n.level > 0 for n in ast.walk(tree)):
            offenders.append(str(path.relative_to(SRC)))
    assert not offenders, f"Imports relatifs interdits (illisibles pour ce test) : {offenders}"


def test_matrix_has_no_cycle() -> None:
    state: dict[str, int] = {}

    def visit(node: str, trail: list[str]) -> None:
        if state.get(node) == 1:
            raise AssertionError("Cycle de dépendances : " + " → ".join([*trail, node]))
        if state.get(node) == 2:
            return
        state[node] = 1
        for dep in ALLOWED.get(node, set()):
            visit(dep, [*trail, node])
        state[node] = 2

    for unit in ALLOWED:
        visit(unit, [])
