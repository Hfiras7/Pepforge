"""Représentation interne du peptide, indépendante de tout champ de force."""

from pepforge.core.peptide import (
    AMBIGUOUS_READING,
    C_TERMINUS_DEFAULTED,
    N_TERMINUS_DEFAULTED,
    Link,
    LinkKind,
    Modification,
    Peptide,
    Residue,
    Terminus,
)
from pepforge.core.structure import Atom, Structure

__all__ = [
    "AMBIGUOUS_READING",
    "C_TERMINUS_DEFAULTED",
    "N_TERMINUS_DEFAULTED",
    "Atom",
    "Link",
    "LinkKind",
    "Modification",
    "Peptide",
    "Residue",
    "Structure",
    "Terminus",
]
