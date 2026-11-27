"""Interface commune aux logiciels de calcul quantique (ORCA, Psi4, Gaussian).

Le projet a besoin de quantique pour dériver et valider les paramètres du
squelette aza-β³ : optimisations de fragments capés, balayages dièdres relaxés
autour de `N–N`, barrière d'inversion de l'azote substitué (sections 15 et 18.8
de `docs/ARCHITECTURE.md`).

Trois décisions portent ce module :

1. **Un seul protocole, plusieurs implémentations.** ORCA et Psi4 sont pris en
   charge derrière `QuantumTool`, ce qui évite d'être prisonnier d'un choix : si
   le laboratoire obtient Gaussian plus tard, c'est une implémentation de plus,
   pas une réécriture.
2. **Le niveau de théorie est une donnée, jamais une valeur par défaut.**
   `QuantumSpec` naît vide et `missing_decisions()` dit ce qu'il reste à fixer.
   Deux entrées de la bibliothèque dérivées à des niveaux différents ne sont pas
   comparables, donc le niveau est fixé une fois pour toute la bibliothèque et
   recopié dans chaque profil produit.
3. **Un profil porte sa provenance.** `TorsionProfile` embarque le logiciel, sa
   version, la spécification de calcul, la date et le matériel. C'est ce qui
   rend un contrôle rejouable des mois plus tard et ce qui sera montré en
   soutenance pour justifier les paramètres.

Module feuille : il ne dépend que de `pepforge.exceptions`, comme le reste de
`pepforge.tools`. Il ne lance aucun calcul et ne lit aucun fichier ; les
implémentations concrètes passent par `ExternalTool` pour que la version soit
enregistrée et la sortie analysée, un calcul non convergé devant lever
`ExternalToolError` plutôt que d'entrer dans une dérivation.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Protocol, runtime_checkable

from pepforge.exceptions import ConfigurationError

HARTREE_TO_KCAL_PER_MOL = 627.509474063
"""Conversion hartree → kcal/mol (constante physique CODATA, pas un paramètre ajusté)."""


class QuantumSoftware(StrEnum):
    """Logiciels pris en charge derrière `QuantumTool`."""

    ORCA = "orca"
    """Recommandé en premier : gratuit pour l'académique sur enregistrement."""

    PSI4 = "psi4"
    """Libre, interface Python donc automatisable et reproductible."""

    GAUSSIAN = "gaussian"
    """Seulement si le laboratoire a déjà la licence."""


@dataclass(frozen=True)
class QuantumSpec:
    """Comment un calcul est fait : logiciel, niveau de théorie, base, balayage.

    Aucune valeur par défaut : le niveau de théorie est une décision
    scientifique à arrêter et à justifier, pas un réglage à deviner.
    """

    software: QuantumSoftware | None = None
    method: str | None = None
    """Méthode, telle qu'elle s'écrit pour le logiciel (ex. `MP2`, `wB97X-D3`)."""

    basis_set: str | None = None
    """Base, telle qu'elle s'écrit pour le logiciel (ex. `6-31G(d)`, `def2-TZVP`)."""

    relaxed_scan: bool | None = None
    """Balayage relaxé (tout optimisé sauf l'angle contraint) ou rigide."""

    scan_step_degrees: float | None = None
    """Pas angulaire du balayage, en degrés."""

    def missing_decisions(self) -> list[str]:
        labels = {
            "software": "logiciel de calcul quantique",
            "method": "méthode",
            "basis_set": "base",
            "relaxed_scan": "balayage relaxé ou rigide",
            "scan_step_degrees": "pas angulaire du balayage",
        }
        return [label for key, label in labels.items() if getattr(self, key) is None]

    def require_complete(self) -> None:
        """Lève `ConfigurationError` tant qu'une décision manque."""
        missing = self.missing_decisions()
        if missing:
            raise ConfigurationError(
                "Spécification de calcul quantique incomplète : " + ", ".join(missing)
            )

    def label(self) -> str:
        """Étiquette courte du niveau de théorie, recopiée dans le rapport."""
        self.require_complete()
        return f"{self.method}/{self.basis_set}"


@dataclass(frozen=True)
class ScanPoint:
    """Un point d'un balayage dièdre : un angle contraint et l'énergie obtenue."""

    angle_degrees: float
    energy_hartree: float


@dataclass(frozen=True)
class TorsionProfile:
    """Profil d'énergie d'une torsion, avec toute sa provenance.

    C'est l'objet de référence auquel le profil du champ de force est comparé
    (contrôle `torsion.profile_vs_quantum`, niveau 3 de la section 18).
    """

    fragment: str
    """Identifiant du fragment capé sur lequel le balayage a été fait."""

    dihedral: tuple[int, int, int, int]
    """Indices des quatre atomes définissant la torsion, dans la géométrie fournie."""

    points: tuple[ScanPoint, ...]
    spec: QuantumSpec
    software_version: str
    """Version exacte du logiciel, telle qu'il l'a déclarée."""

    computed_on: str = ""
    """Date du calcul, au format ISO."""

    hardware: str = ""
    """Matériel utilisé, pour que deux profils restent comparables."""

    metadata: dict[str, str] = field(default_factory=dict)

    def relative_energies_kcal(self) -> tuple[float, ...]:
        """Énergies relatives au minimum du profil, en kcal/mol.

        Le profil n'a de sens qu'en relatif : seules les différences d'énergie
        se comparent à celles du champ de force.
        """
        if not self.points:
            return ()
        floor = min(point.energy_hartree for point in self.points)
        return tuple(
            (point.energy_hartree - floor) * HARTREE_TO_KCAL_PER_MOL for point in self.points
        )


@runtime_checkable
class QuantumTool(Protocol):
    """Un logiciel de calcul quantique, vu par PepForge.

    Toute implémentation passe sa fonction de contrat
    `tests/contracts.py::assert_quantum_tool_contract`.
    """

    @property
    def software(self) -> QuantumSoftware:
        """Logiciel implémenté."""
        ...

    def version(self) -> str:
        """Version déclarée par l'exécutable, enregistrée dans le rapport."""
        ...

    def optimize(self, geometry: Path, spec: QuantumSpec, outdir: Path) -> Path:
        """Optimise une géométrie et retourne le chemin de la géométrie obtenue."""
        ...

    def scan_dihedral(
        self,
        geometry: Path,
        dihedral: Sequence[int],
        spec: QuantumSpec,
        outdir: Path,
    ) -> TorsionProfile:
        """Balaie une torsion sur un tour complet et retourne le profil obtenu."""
        ...

    def charges(self, geometry: Path, spec: QuantumSpec, outdir: Path) -> dict[str, float]:
        """Charges atomiques, indexées par nom d'atome."""
        ...
