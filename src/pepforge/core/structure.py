"""Structure 3D du soluté, indépendante de tout champ de force et de tout moteur.

`Structure` est le type échangé entre `io` (lecture PDB/mmCIF), `build`
(construction et correction 3D) et `engines` (écriture de la topologie). Les
noms d'atomes et de résidus y suivent la nomenclature du format d'origine
(PDB) ; la traduction vers la nomenclature d'un champ de force (résidus
terminaux, tautomères, patches) relève exclusivement de `pepforge.forcefields`.

Les coordonnées sont en angströms, comme dans les formats PDB et mmCIF.

Identification des résidus (étape 3) : un résidu est désigné par la chaîne, le
numéro et le code d'insertion **de l'auteur** (`auth_asym_id`, `auth_seq_id`,
`pdbx_PDB_ins_code` en mmCIF ; colonnes 22, 23-26 et 27 en PDB), c'est-à-dire
la numérotation que montrent les articles et les fichiers PDB. Deux résidus de
même numéro et de codes d'insertion différents (`52` et `52A`) sont distincts.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from pepforge.core.peptide import Peptide

STRUCTURE_FORMAT = "structure_format"
"""Clé de `Structure.metadata` : format du fichier lu (`pdb`, `mmcif`)."""

STRUCTURE_SHA256 = "structure_sha256"
"""Clé de `Structure.metadata` : empreinte SHA-256 du fichier lu, tel qu'il a été reçu."""

STRUCTURE_MODEL_COUNT = "structure_model_count"
"""Clé de `Structure.metadata` : nombre de modèles présents dans le fichier lu."""

STRUCTURE_DECLARED_MODEL_COUNT = "structure_declared_model_count"
"""Clé de `Structure.metadata` : nombre de modèles annoncé par l'en-tête (`NUMMDL`), s'il existe."""


@dataclass(frozen=True, order=True)
class ResidueId:
    """Identifiant d'un résidu dans une structure : chaîne, numéro, code d'insertion, nom."""

    chain_id: str
    number: int
    insertion_code: str = ""
    """Code d'insertion (`A` dans `52A`) ; chaîne vide s'il n'y en a pas."""

    name: str = ""

    def label(self) -> str:
        """Forme lisible : `A:ARG52A`."""
        return f"{self.chain_id}:{self.name}{self.number}{self.insertion_code}"


@dataclass(frozen=True)
class MissingAtom:
    """Atome que le fichier déclare non observé (`REMARK 470`, `pdbx_unobs_or_zero_occ_atoms`)."""

    residue: ResidueId
    atom: str

    def label(self) -> str:
        return f"{self.residue.label()}:{self.atom}"


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

    insertion_code: str = ""
    """Code d'insertion du résidu ; chaîne vide s'il n'y en a pas."""

    occupancy: float | None = None
    """Occupation lue dans le fichier ; `None` pour un atome qui n'en vient pas (construit)."""

    b_factor: float | None = None
    """Facteur de température lu dans le fichier ; `None` s'il n'en vient pas."""

    def __post_init__(self) -> None:
        for axis, value in (("x", self.x), ("y", self.y), ("z", self.z)):
            if not math.isfinite(value):
                raise ValueError(f"Atome {self.serial} ({self.name}) : coordonnée {axis} non finie")

    @property
    def residue_id(self) -> ResidueId:
        return ResidueId(self.chain_id, self.residue_number, self.insertion_code, self.residue_name)


class ConnectionKind:
    """Natures de liaison déclarées dans un fichier (`SSBOND`, `LINK`, `struct_conn`)."""

    DISULFIDE = "disulfide"
    COVALENT = "covalent"
    METAL = "metal_coordination"
    HYDROGEN = "hydrogen_bond"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class Connection:
    """Liaison déclarée par le fichier entre deux atomes (pont disulfure, liaison covalente...).

    Ce n'est pas une liaison déduite des distances : c'est ce que le fichier
    annonce, gardé pour que rien ne disparaisse en silence à la lecture.
    """

    kind: str
    first: ResidueId
    first_atom: str
    second: ResidueId
    second_atom: str


def _empty_metadata() -> Mapping[str, str]:
    return MappingProxyType({})


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

    metadata: Mapping[str, str] = field(default_factory=_empty_metadata)

    model: int | None = None
    """Numéro du modèle lu (fichier multi-modèles, RMN) ; `None` si la structure n'en vient pas."""

    connections: tuple[Connection, ...] = ()
    """Liaisons déclarées par le fichier (`SSBOND`, `LINK`, `struct_conn`)."""

    declared_missing_residues: tuple[ResidueId, ...] = ()
    """Résidus déclarés non observés (`REMARK 465`, `pdbx_unobs_or_zero_occ_residues`)."""

    declared_missing_atoms: tuple[MissingAtom, ...] = ()
    """Atomes déclarés non observés dans des résidus présents (`REMARK 470`, mmCIF)."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def __len__(self) -> int:
        return len(self.atoms)

    @property
    def chain_ids(self) -> list[str]:
        """Identifiants de chaînes, dans l'ordre d'apparition."""
        return list(dict.fromkeys(atom.chain_id for atom in self.atoms))

    def residue_ids(self) -> list[ResidueId]:
        """Résidus, dans l'ordre d'apparition des atomes."""
        return list(dict.fromkeys(atom.residue_id for atom in self.atoms))

    def has_alternate_locations(self) -> bool:
        return any(atom.alt_loc for atom in self.atoms)

    def has_hydrogens(self) -> bool:
        return any(atom.element.upper() in ("H", "D") for atom in self.atoms)
