"""Représentation interne du peptide, indépendante de tout champ de force."""

from pepforge.core.peptide import Peptide, Residue, Terminus
from pepforge.core.structure import Atom, Structure

__all__ = ["Atom", "Peptide", "Residue", "Structure", "Terminus"]
