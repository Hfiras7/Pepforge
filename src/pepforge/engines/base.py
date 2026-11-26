"""Interface commune des générateurs de fichiers moteur."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from pepforge.config import Engine, PrepConfig
from pepforge.forcefields.base import ParameterSet


@runtime_checkable
class EngineWriter(Protocol):
    """Contrat que respecte chaque générateur de fichiers d'entrée."""

    engine: Engine

    def write(self, parameters: ParameterSet, config: PrepConfig, outdir: Path) -> list[Path]:
        """Écrit les fichiers d'entrée et retourne leurs chemins.

        Toute conversion depuis un autre moteur est une étape séparée et vérifiée,
        jamais réalisée implicitement ici.
        """
        ...
