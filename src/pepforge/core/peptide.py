"""Modèle de données du peptide.

Cette représentation ne contient aucun nom de résidu propre à un champ de
force (par exemple `NHE` ou `NME` d'AMBER) : la traduction se fait dans
`pepforge.forcefields`.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Terminus:
    """Extrémité N- ou C-terminale (libre, acétylée, amidée...)."""

    key: str
    """Identifiant interne stable, par exemple `free_amine` ou `amide`."""

    label: str
    """Libellé lisible en français."""


@dataclass(frozen=True)
class Residue:
    """Résidu de la chaîne peptidique."""

    code: str
    """Code à trois lettres (ou code de la bibliothèque pour un non standard)."""

    one_letter: str | None = None
    standard: bool = True


@dataclass
class Peptide:
    """Peptide linéaire : chaîne de résidus et deux extrémités explicites."""

    residues: list[Residue]
    n_terminus: Terminus
    c_terminus: Terminus
    metadata: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.residues)

    @property
    def one_letter_sequence(self) -> str:
        """Séquence à une lettre ; `X` pour un résidu sans code à une lettre."""
        return "".join(r.one_letter or "X" for r in self.residues)

    @property
    def three_letter_sequence(self) -> str:
        return "-".join(r.code for r in self.residues)

    def composition(self) -> dict[str, int]:
        """Nombre d'occurrences de chaque résidu (code à trois lettres)."""
        return dict(Counter(r.code for r in self.residues))
