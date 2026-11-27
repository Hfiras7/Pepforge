"""Interface commune des adaptateurs de champ de force.

Un adaptateur traduit le peptide de la représentation interne (`pepforge.core`)
vers la nomenclature et les jeux de paramètres d'un champ de force. Il
n'écrit jamais de fichier moteur : il produit un `ParameterSet`, que
`pepforge.engines` sérialise.

Contrat (vérifié par `tests/test_interfaces.py` et, pour chaque adaptateur
réel, par un test de charge nette et de paramètres manquants) :

- entrée : un `Peptide` dont les extrémités sont explicites ;
- sortie : un `ParameterSet` dont `missing` est vide ;
- un paramètre absent lève `MissingParameterError` : il n'est jamais complété
  par une valeur inventée ;
- l'adaptateur ne connaît pas le moteur, sauf pour déclarer ceux qu'il prend
  en charge (`supported_engines`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, runtime_checkable

from pepforge.config import Engine
from pepforge.core.peptide import Peptide


@dataclass(frozen=True)
class ParameterSource:
    """Provenance d'un jeu de paramètres."""

    name: str
    """Nom du jeu de paramètres (ex. un fichier `frcmod`, un `stream file` CHARMM)."""

    citation: str
    """Publication ou procédure de dérivation documentée."""

    path: Path | None = None
    """Fichier de paramètres à charger par le moteur, s'il y en a un."""

    version: str | None = None


@dataclass
class ParameterSet:
    """Résultat de la paramétrisation d'un peptide pour un champ de force donné."""

    forcefield: str
    residue_names: list[str]
    """Noms des résidus dans la nomenclature du champ de force (ex. `NHE` en AMBER),
    dans l'ordre de la chaîne, extrémités comprises lorsqu'elles sont des résidus."""

    sources: dict[str, str] = field(default_factory=dict)
    """Provenance de chaque jeu de paramètres (publication, fichier, version)."""

    parameter_files: list[ParameterSource] = field(default_factory=list)
    """Fichiers supplémentaires à charger (résidus non standards, fragments)."""

    net_charge: float | None = None
    """Charge nette d'après les charges partielles du champ de force."""

    missing: list[str] = field(default_factory=list)
    """Paramètres absents : doit être vide avant toute écriture de fichier moteur."""


@runtime_checkable
class ForceFieldAdapter(Protocol):
    """Contrat que respecte chaque adaptateur de champ de force."""

    name: str
    supported_engines: frozenset[Engine]
    compatible_water_models: frozenset[str]
    """Modèles d'eau publiés comme compatibles avec ce champ de force."""

    def parametrize(self, peptide: Peptide) -> ParameterSet:
        """Associe les paramètres du champ de force au peptide.

        Doit lever `MissingParameterError` plutôt que de compléter un paramètre
        absent par une valeur inventée.
        """
        ...
