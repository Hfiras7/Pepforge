"""Encapsulation des exécutables externes (tleap, antechamber, parmchk2, gmx, psfgen).

La version de l'outil est enregistrée, et la sortie est analysée à la recherche
de motifs d'échec : le code de retour seul ne suffit pas (tleap, par exemple,
retourne 0 même quand des paramètres manquent).
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from pepforge.exceptions import ExternalToolError


@dataclass
class ToolResult:
    """Résultat d'un appel à un outil externe."""

    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    version: str | None


@dataclass
class ExternalTool:
    """Un exécutable externe et les motifs qui signalent un échec dans sa sortie."""

    executable: str
    version_args: Sequence[str] = ("--version",)
    failure_patterns: Sequence[str] = field(default_factory=tuple)

    def path(self) -> str:
        found = shutil.which(self.executable)
        if found is None:
            raise ExternalToolError(f"Exécutable introuvable dans le PATH : {self.executable}")
        return found

    def version(self) -> str | None:
        """Première ligne non vide de la sortie de la commande de version."""
        proc = subprocess.run(
            [self.path(), *self.version_args], capture_output=True, text=True, check=False
        )
        for line in (proc.stdout + proc.stderr).splitlines():
            if line.strip():
                return line.strip()
        return None

    def run(self, args: Sequence[str], cwd: Path | None = None) -> ToolResult:
        """Exécute l'outil, puis vérifie le code de retour et la sortie."""
        command = [self.path(), *args]
        proc = subprocess.run(command, capture_output=True, text=True, cwd=cwd, check=False)
        result = ToolResult(
            command=command,
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            version=self.version(),
        )
        if proc.returncode != 0:
            raise ExternalToolError(
                f"{self.executable} a échoué (code {proc.returncode}) : {proc.stderr.strip()}"
            )
        output = proc.stdout + proc.stderr
        for pattern in self.failure_patterns:
            match = re.search(pattern, output, flags=re.MULTILINE)
            if match:
                raise ExternalToolError(
                    f"{self.executable} : motif d'échec détecté dans la sortie : {match.group(0)!r}"
                )
        return result
