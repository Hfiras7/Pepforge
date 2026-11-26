"""Encapsulation des exécutables externes (tleap, antechamber, parmchk2, gmx, psfgen).

Module feuille : il ne dépend que de `pepforge.exceptions`, ce qui permet aux
adaptateurs de champ de force, aux constructeurs 3D et aux générateurs moteur
de l'utiliser sans créer de dépendance circulaire avec `pepforge.pipeline`.
"""

from pepforge.tools.external import ExternalTool, ToolResult

__all__ = ["ExternalTool", "ToolResult"]
