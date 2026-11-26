"""Chargement de la bibliothèque de résidus depuis les fichiers YAML de `data/`.

Chaque fichier de données déclare sa source : aucune donnée n'est codée en dur
dans la logique.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from typing import Any

import yaml

from pepforge.core.peptide import Residue, Terminus
from pepforge.exceptions import UnknownResidueError, UnknownTerminusError


@dataclass
class ResidueLibrary:
    """Résidus et extrémités disponibles, avec la source de chaque fichier chargé."""

    residues: dict[str, Residue]
    n_termini: dict[str, Terminus]
    c_termini: dict[str, Terminus]
    n_aliases: dict[str, str]
    c_aliases: dict[str, str]
    default_n: str
    default_c: str
    sources: dict[str, str]

    def by_code(self, code: str) -> Residue:
        try:
            return self.residues[code.upper()]
        except KeyError:
            raise UnknownResidueError(f"Résidu inconnu : {code!r}") from None

    def by_one_letter(self, letter: str) -> Residue:
        for residue in self.residues.values():
            if residue.one_letter == letter.upper():
                return residue
        raise UnknownResidueError(f"Code à une lettre inconnu : {letter!r}")

    def has_code(self, code: str) -> bool:
        return code.upper() in self.residues

    def n_terminus(self, alias: str | None = None) -> Terminus:
        if alias is None:
            return self.n_termini[self.default_n]
        key = self.n_aliases.get(alias.lower())
        if key is None:
            raise UnknownTerminusError(f"Extrémité N-terminale inconnue : {alias!r}")
        return self.n_termini[key]

    def c_terminus(self, alias: str | None = None) -> Terminus:
        if alias is None:
            return self.c_termini[self.default_c]
        key = self.c_aliases.get(alias.lower())
        if key is None:
            raise UnknownTerminusError(f"Extrémité C-terminale inconnue : {alias!r}")
        return self.c_termini[key]

    def is_n_alias(self, token: str) -> bool:
        return token.lower() in self.n_aliases

    def is_c_alias(self, token: str) -> bool:
        return token.lower() in self.c_aliases


def _read_yaml(name: str) -> dict[str, Any]:
    text = resources.files("pepforge.residues").joinpath("data", name).read_text("utf-8")
    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"Fichier de données invalide : {name}")
    return data


def _termini(entries: list[dict[str, Any]]) -> tuple[dict[str, Terminus], dict[str, str]]:
    termini: dict[str, Terminus] = {}
    aliases: dict[str, str] = {}
    for entry in entries:
        terminus = Terminus(key=entry["key"], label=entry["label"])
        termini[terminus.key] = terminus
        for alias in entry["aliases"]:
            aliases[str(alias).lower()] = terminus.key
    return termini, aliases


def load_default_library() -> ResidueLibrary:
    """Charge les acides aminés standards et les extrémités usuelles."""
    standard = _read_yaml("standard_amino_acids.yaml")
    termini = _read_yaml("termini.yaml")
    residues = {
        entry["code"]: Residue(code=entry["code"], one_letter=entry["one_letter"], standard=True)
        for entry in standard["residues"]
    }
    n_termini, n_aliases = _termini(termini["n_terminus"])
    c_termini, c_aliases = _termini(termini["c_terminus"])
    return ResidueLibrary(
        residues=residues,
        n_termini=n_termini,
        c_termini=c_termini,
        n_aliases=n_aliases,
        c_aliases=c_aliases,
        default_n=termini["default_n_terminus"],
        default_c=termini["default_c_terminus"],
        sources={
            "standard_amino_acids.yaml": standard["source"],
            "termini.yaml": termini["source"],
        },
    )
