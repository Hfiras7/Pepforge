"""Interface commune des constructeurs et correcteurs de structure 3D.

Deux familles de composants, choisis indépendamment du champ de force :

- `StructureBuilder` : construit une structure à partir d'un `Peptide`
  (conformation étendue ou hélicoïdale, hydrogènes compris) ;
- `StructureModifier` : transforme une structure existante (mutation,
  changement d'extrémité, reconstruction d'atomes manquants, protonation).

Contrat commun :

- la structure d'entrée n'est jamais modifiée sur place (`Structure` est
  immuable) : une nouvelle structure est retournée ;
- chaque transformation est enregistrée dans le `Journal` reçu : quoi,
  pourquoi, avec quels paramètres ;
- les extrémités du `Peptide` sont conservées telles quelles ;
- une incohérence lève `StructureError` ou `ChemistryError`.

L'outil de construction (tleap, un autre constructeur...) n'est pas encore
choisi : il sera encapsulé par `pepforge.tools.ExternalTool` si c'est un
exécutable externe.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from pepforge.core.peptide import Peptide
from pepforge.core.structure import Structure
from pepforge.report.journal import Journal


@runtime_checkable
class StructureBuilder(Protocol):
    """Construit une structure 3D initiale à partir d'un peptide."""

    name: str

    def build(self, peptide: Peptide, journal: Journal) -> Structure:
        """Retourne une structure dont `peptide` est l'argument reçu."""
        ...


@runtime_checkable
class StructureModifier(Protocol):
    """Transforme une structure existante et journalise la transformation."""

    name: str

    def apply(self, structure: Structure, journal: Journal) -> Structure:
        """Retourne une nouvelle structure ; `structure` reste inchangée."""
        ...
