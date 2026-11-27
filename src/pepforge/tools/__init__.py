"""Encapsulation des outils externes : exécutables classiques et calcul quantique.

Module feuille : il ne dépend que de `pepforge.exceptions`, ce qui permet aux
adaptateurs de champ de force, aux constructeurs 3D et aux générateurs moteur
de l'utiliser sans créer de dépendance circulaire avec `pepforge.pipeline`.
"""

from pepforge.tools.external import ExternalTool, ToolResult
from pepforge.tools.quantum import (
    HARTREE_TO_KCAL_PER_MOL,
    QuantumSoftware,
    QuantumSpec,
    QuantumTool,
    ScanPoint,
    TorsionProfile,
)

__all__ = [
    "HARTREE_TO_KCAL_PER_MOL",
    "ExternalTool",
    "QuantumSoftware",
    "QuantumSpec",
    "QuantumTool",
    "ScanPoint",
    "ToolResult",
    "TorsionProfile",
]
