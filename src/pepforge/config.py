"""Configuration d'une préparation.

Le moteur, le champ de force, le modèle d'eau et l'objectif de simulation ne
sont pas encore arrêtés pour le projet : ils restent optionnels ici, et aucune
valeur n'est supposée par défaut. `missing_decisions()` liste ce qui manque
avant de pouvoir lancer une préparation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import IntEnum, StrEnum
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


class ValidationLevel(IntEnum):
    """Jusqu'où pousser la vérification d'une sortie (`docs/ARCHITECTURE.md`, section 18).

    Les niveaux sont **cumulatifs** : demander le niveau 2 exige aussi les
    niveaux 0 et 1. Le coût croît d'un facteur considérable à chaque palier,
    d'où le choix laissé à l'utilisateur. Un niveau demandé dont un contrôle
    n'a pas pu être exécuté est rapporté en échec, jamais passé sous silence
    (`pepforge.checks.levels.audit`).
    """

    FILE_CONFORMITY = 0
    """Conformité des fichiers : le fichier dit-il ce que la représentation interne disait ?

    Gratuit, automatisable, entièrement hors ligne. C'est ce palier qui attrape
    les erreurs silencieuses du catalogue de la section 17.2.
    """

    NUMERIC_COHERENCE = 1
    """Cohérence numérique : énergie en point unique, minimisation courte qui converge."""

    STRUCTURAL_COHERENCE = 2
    """Cohérence structurale : dynamique courte, géométrie comparée au fragment quantique."""

    ENERGETIC_COHERENCE = 3
    """Cohérence énergétique : profils de torsion contre quantique, répliques, peptide témoin."""

    EXPERIMENTAL_AGREEMENT = 4
    """Accord expérimental (RMN, affinité). Hors de l'outil : la confrontation est manuelle."""


@dataclass(frozen=True)
class SolvationSpec:
    """Paramètres de solvatation et de neutralisation, tous à décider.

    Aucune valeur par défaut : la marge de boîte et la concentration en sel
    dépendent de l'objectif (peptide libre ou complexe) et doivent être
    justifiées dans la configuration (`PrepConfig.notes`).
    """

    box_shape: str | None = None
    """Forme de la boîte périodique (ex. octaèdre tronqué, cube)."""

    margin_angstrom: float | None = None
    """Distance minimale soluté–bord de boîte, en angströms."""

    salt_concentration_molar: float | None = None
    """Concentration en sel visée, en mol/L."""

    neutralize: bool | None = None

    def missing_decisions(self) -> list[str]:
        labels = {
            "box_shape": "forme de la boîte",
            "margin_angstrom": "marge de la boîte",
            "salt_concentration_molar": "concentration en sel",
            "neutralize": "neutralisation",
        }
        return [label for key, label in labels.items() if getattr(self, key) is None]


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
    protonation: dict[str, str] = field(default_factory=dict)
    """États de protonation imposés, qui l'emportent sur le choix au pH.

    Clé : `"N"` ou `"C"` pour une extrémité, ou la position du résidu (à partir
    de 1) ; valeur : clé interne de l'état (`his_epsilon`...). Vide par défaut,
    ce qui ne suppose rien : chaque état non imposé est déduit du pH et de son
    pKa sourcé. Obligatoire pour une histidine neutre, dont le pKa ne départage
    pas les deux tautomères. La raison de chaque choix s'écrit dans `notes`.
    """

    receptor: str | None = None
    solvation: SolvationSpec = field(default_factory=SolvationSpec)
    """Vérifiée par l'étape de solvatation elle-même (`SolvationSpec.missing_decisions`)."""

    validation_level: ValidationLevel = ValidationLevel.FILE_CONFORMITY
    """Seul champ avec une valeur par défaut, parce que le niveau 0 ne coûte rien.

    Les paliers supérieurs demandent des simulations ou du calcul quantique :
    ils sont choisis explicitement, et l'outil refuse de les déclarer satisfaits
    tant qu'il ne les a pas réellement exécutés.
    """

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
        data["validation_level"] = int(self.validation_level)
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
            if values.get("validation_level") is not None:
                values["validation_level"] = ValidationLevel(values["validation_level"])
            if values.get("protonation") is None:
                values["protonation"] = {}
            if not isinstance(values["protonation"], dict):
                raise ConfigurationError("`protonation` doit être un dictionnaire")
            values["protonation"] = {str(k): v for k, v in values["protonation"].items()}
            if "solvation" in values:
                solvation = values["solvation"] or {}
                if not isinstance(solvation, dict):
                    raise ConfigurationError("`solvation` doit être un dictionnaire")
                values["solvation"] = SolvationSpec(**solvation)
        except (ValueError, TypeError) as exc:
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
