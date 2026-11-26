"""Adaptateurs de champ de force (AMBER, CHARMM, GROMOS).

Aucun champ de force n'est encore choisi : seule l'interface commune est définie.
"""

from pepforge.forcefields.base import ForceFieldAdapter, ParameterSet

__all__ = ["ForceFieldAdapter", "ParameterSet"]
