"""Interface commune des générateurs de fichiers d'entrée par moteur.

Un `EngineWriter` sérialise, pour un moteur donné, une `Structure` et le
`ParameterSet` produit par un adaptateur de champ de force. C'est le seul
composant qui écrit des fichiers propres à un moteur (`.prmtop`, `.top`,
`.psf`...).

Contrat :

- refuser un `ParameterSet` dont `missing` n'est pas vide
  (`MissingParameterError`) ;
- refuser un champ de force qui ne déclare pas ce moteur dans
  `supported_engines` (`ConfigurationError`) ;
- ne jamais convertir implicitement depuis un autre moteur : une conversion
  est une étape séparée du pipeline, suivie de ses propres contrôles ;
- passer tout appel à `tleap`, `gmx`, `psfgen`... par
  `pepforge.tools.ExternalTool`, et reporter les versions dans
  `EngineFiles.tool_versions`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, runtime_checkable

from pepforge.config import Engine, SolvationSpec
from pepforge.core.structure import Structure
from pepforge.forcefields.base import ParameterSet
from pepforge.report.journal import Journal


@dataclass(frozen=True)
class EngineFiles:
    """Fichiers produits pour un moteur, indexés par leur rôle."""

    engine: Engine
    files: dict[str, Path]
    """Rôle → chemin, par exemple `topology`, `coordinates`, `parameters`."""

    solvated: bool = False
    net_charge: float | None = None
    """Charge nette relue dans les fichiers produits (et non recopiée de l'entrée)."""

    atom_count: int | None = None
    tool_versions: dict[str, str] = field(default_factory=dict)


@runtime_checkable
class EngineWriter(Protocol):
    """Contrat que respecte chaque générateur de fichiers d'entrée."""

    engine: Engine

    def write_topology(
        self, structure: Structure, parameters: ParameterSet, outdir: Path
    ) -> EngineFiles:
        """Écrit topologie, paramètres et coordonnées du soluté dans le vide."""
        ...

    def solvate(
        self, files: EngineFiles, spec: SolvationSpec, outdir: Path, journal: Journal
    ) -> EngineFiles:
        """Solvate, neutralise et ajoute le sel selon `spec`, puis journalise."""
        ...
