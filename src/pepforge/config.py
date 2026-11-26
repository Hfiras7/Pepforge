"""Configuration d'une préparation.

Le moteur, le champ de force, le modèle d'eau et l'objectif de simulation ne
sont pas encore arrêtés pour le projet : ils restent optionnels ici, et aucune
valeur n'est supposée par défaut. `missing_decisions()` liste ce qui manque
avant de pouvoir lancer une préparation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any

import yaml

from pepforge.exceptions import ConfigurationError


class Engine(StrEnum):
    """Moteurs de dynamique moléculaire visés."""

    AMBER = "amber"
    GROMACS = "gromacs"
    NAMD = "namd"


class SimulationGoal(StrEnum):
    """Les deux scénarios de simulation, à ne jamais confondre."""

    FREE_PEPTIDE = "free_peptide"
    """Peptide libre en solution : stabilité conformationnelle."""

    PEPTIDE_RECEPTOR_COMPLEX = "peptide_receptor_complex"
    """Complexe peptide–récepteur (DCC, TrkC, UNC5)."""


@dataclass
class PrepConfig:
    """Paramètres d'une préparation, sérialisables en YAML pour la reproductibilité."""

    sequence: str | None = None
    input_structure: str | None = None
    engine: Engine | None = None
    forcefield: str | None = None
    water_model: str | None = None
    goal: SimulationGoal | None = None
    ph: float | None = None
    receptor: str | None = None
    notes: dict[str, str] = field(default_factory=dict)

    def missing_decisions(self) -> list[str]:
        """Retourne la liste des choix encore à faire avant une préparation."""
        missing: list[str] = []
        if self.sequence is None and self.input_structure is None:
            missing.append("entrée (séquence ou structure 3D)")
        if self.engine is None:
            missing.append("moteur de simulation")
        if self.forcefield is None:
            missing.append("champ de force")
        if self.water_model is None:
            missing.append("modèle d'eau")
        if self.goal is None:
            missing.append("objectif (peptide libre ou complexe peptide–récepteur)")
        if self.ph is None:
            missing.append("pH de référence")
        if self.goal is SimulationGoal.PEPTIDE_RECEPTOR_COMPLEX and self.receptor is None:
            missing.append("récepteur (structure ou code PDB)")
        return missing

    def require_complete(self) -> None:
        """Lève `ConfigurationError` si un choix manque."""
        missing = self.missing_decisions()
        if missing:
            raise ConfigurationError("Configuration incomplète : " + ", ".join(missing))

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["engine"] = self.engine.value if self.engine else None
        data["goal"] = self.goal.value if self.goal else None
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PrepConfig:
        known = set(cls.__dataclass_fields__)
        unknown = set(data) - known
        if unknown:
            raise ConfigurationError(f"Clés inconnues dans la configuration : {sorted(unknown)}")
        values = dict(data)
        try:
            if values.get("engine") is not None:
                values["engine"] = Engine(values["engine"])
            if values.get("goal") is not None:
                values["goal"] = SimulationGoal(values["goal"])
        except ValueError as exc:
            raise ConfigurationError(str(exc)) from exc
        return cls(**values)

    def to_yaml(self, path: Path) -> None:
        path.write_text(
            yaml.safe_dump(self.to_dict(), allow_unicode=True, sort_keys=False), encoding="utf-8"
        )

    @classmethod
    def from_yaml(cls, path: Path) -> PrepConfig:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        if not isinstance(data, dict):
            raise ConfigurationError(f"{path} ne contient pas un dictionnaire YAML")
        return cls.from_dict(data)
