"""Composants chimiques reconnus à la lecture d'une structure (eau, groupes terminaux).

Les données viennent de `data/structure_components.yaml`, avec leur source
(dictionnaire des composants chimiques de la wwPDB). Le fichier est vérifié au
chargement : chaque groupe terminal renvoie à une extrémité connue de
`termini.yaml`, sinon `ReferenceDataError`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pepforge.exceptions import ReferenceDataError
from pepforge.residues.library import ResidueLibrary, load_default_library, read_data_file

DATA_FILE = "structure_components.yaml"


@dataclass(frozen=True)
class StructureComponents:
    """Noms de composants que l'inventaire et la correspondance savent classer."""

    water: frozenset[str]
    n_terminal_caps: dict[str, str]
    """Composant (`ACE`) → clé d'extrémité N (`acetyl`)."""

    c_terminal_caps: dict[str, str]
    """Composant (`NH2`, `NME`) → clé d'extrémité C (`amide`, `n_methylamide`)."""

    free_c_terminus: str
    """Extrémité C prouvée par un atome de `free_c_terminus_atoms` (carboxyle libre)."""

    free_c_terminus_atoms: frozenset[str]
    free_n_terminus: str
    """Extrémité N prouvée par un hydrogène de `free_n_terminus_hydrogens` (amine libre)."""

    free_n_terminus_hydrogens: frozenset[str]
    source: str


def _names(data: dict[str, Any], key: str) -> frozenset[str]:
    value = data.get(key)
    if not isinstance(value, list) or not value:
        raise ReferenceDataError(f"{DATA_FILE} : liste « {key} » absente ou vide")
    return frozenset(str(v).upper() for v in value)


def _caps(data: dict[str, Any], key: str, known: dict[str, Any]) -> dict[str, str]:
    caps: dict[str, str] = {}
    for entry in data.get(key) or []:
        component = str(entry["component"]).upper()
        terminus = str(entry["terminus"])
        if terminus not in known:
            raise ReferenceDataError(
                f"{DATA_FILE} : le groupe {component} renvoie à l'extrémité inconnue {terminus!r}"
            )
        caps[component] = terminus
    return caps


def _terminus(data: dict[str, Any], key: str, known: dict[str, Any]) -> str:
    value = data.get(key)
    if value not in known:
        raise ReferenceDataError(
            f"{DATA_FILE} : « {key} » renvoie à l'extrémité inconnue {value!r}"
        )
    return str(value)


def load_structure_components(library: ResidueLibrary | None = None) -> StructureComponents:
    """Charge et vérifie `structure_components.yaml`."""
    library = library or load_default_library()
    data = read_data_file(DATA_FILE)
    source = data.get("source")
    if not source:
        raise ReferenceDataError(f"{DATA_FILE} : source non déclarée")
    return StructureComponents(
        water=_names(data, "water"),
        n_terminal_caps=_caps(data, "n_terminal_caps", dict(library.n_termini)),
        c_terminal_caps=_caps(data, "c_terminal_caps", dict(library.c_termini)),
        free_c_terminus=_terminus(data, "free_c_terminus", dict(library.c_termini)),
        free_c_terminus_atoms=_names(data, "free_c_terminus_atoms"),
        free_n_terminus=_terminus(data, "free_n_terminus", dict(library.n_termini)),
        free_n_terminus_hydrogens=_names(data, "free_n_terminus_hydrogens"),
        source=str(source),
    )
