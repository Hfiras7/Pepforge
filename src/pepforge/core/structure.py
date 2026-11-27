"""Structure 3D du soluté, indépendante de tout champ de force et de tout moteur.

`Structure` est le type échangé entre `io` (lecture PDB/mmCIF), `build`
(construction et correction 3D) et `engines` (écriture de la topologie). Les
noms d'atomes et de résidus y suivent la nomenclature du format d'origine
(PDB) ; la traduction vers la nomenclature d'un champ de force (résidus
terminaux, tautomères, patches) relève exclusivement de `pepforge.forcefields`.

Les coordonnées sont en angströms, comme dans les formats PDB et mmCIF.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pepforge.core.peptide import Peptide


@dataclass(frozen=True)
class Atom:
    """Atome d'une structure 3D."""

    serial: int
    name: str
    element: str
    residue_name: str
    residue_number: int
    chain_id: str
    x: float
    y: float
    z: float
    alt_loc: str | None = None
    """Indicateur de conformation alternative, à résoudre avant toute paramétrisation."""

    hetero: bool = False
    """Vrai pour un enregistrement HETATM (ligand, ion, eau, résidu non standard)."""


@dataclass(frozen=True)
class Structure:
    """Ensemble d'atomes, éventuellement associé au peptide qu'il représente.

    Immuable : une opération qui modifie la structure en produit une nouvelle
    et enregistre l'opération dans le journal (`pepforge.report.Journal`).
    """

    atoms: tuple[Atom, ...]
    peptide: Peptide | None = None
    """Peptide décrit par cette structure, si la correspondance a été établie."""

    source: str = ""
    """Origine : fichier lu, constructeur utilisé, ou étape qui l'a produite."""

    metadata: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.atoms)

    @property
    def chain_ids(self) -> list[str]:
        """Identifiants de chaînes, dans l'ordre d'apparition."""
        return list(dict.fromkeys(atom.chain_id for atom in self.atoms))

    def has_alternate_locations(self) -> bool:
        return any(atom.alt_loc for atom in self.atoms)

    def has_hydrogens(self) -> bool:
        return any(atom.element.upper() == "H" for atom in self.atoms)
