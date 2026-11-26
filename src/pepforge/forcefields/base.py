"""Interface commune des adaptateurs de champ de force.

Un adaptateur n'écrit jamais de fichier moteur : il produit un `ParameterSet`,
que `pepforge.engines` sérialise.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

from pepforge.config import Engine
from pepforge.core.peptide import Peptide


@dataclass
class ParameterSet:
    """Résultat de la paramétrisation d'un peptide pour un champ de force donné."""

    forcefield: str
    residue_names: list[str]
    """Noms des résidus dans la nomenclature du champ de force (ex. `NHE` en AMBER)."""

    sources: dict[str, str] = field(default_factory=dict)
    """Provenance de chaque jeu de paramètres (publication, fichier, version)."""

    missing: list[str] = field(default_factory=list)
    """Paramètres absents : doit être vide avant toute écriture de fichier moteur."""


@runtime_checkable
class ForceFieldAdapter(Protocol):
    """Contrat que respecte chaque adaptateur de champ de force."""

    name: str
    supported_engines: frozenset[Engine]

    def parametrize(self, peptide: Peptide) -> ParameterSet:
        """Associe les paramètres du champ de force au peptide.

        Doit lever `MissingParameterError` plutôt que de compléter un paramètre
        absent par une valeur inventée.
        """
        ...
