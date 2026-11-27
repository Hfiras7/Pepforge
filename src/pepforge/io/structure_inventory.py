"""Inventaire d'une structure importée, et correspondance `Structure` → `Peptide`.

L'inventaire décrit sans rien modifier : chaînes, résidus, eaux, ions,
composants non reconnus, hydrogènes, conformations alternatives, codes
d'insertion, trous de numérotation, résidus déclarés non observés, liaisons
déclarées. Il se termine par la liste des **décisions à prendre** avant
d'utiliser la structure (`pending_decisions`) : rien n'est tranché à la place
de l'utilisateur.

Règles de classement, toutes structurelles ou tirées d'un fichier de données
sourcé, jamais d'une liste écrite dans le code :

- **acide aminé standard** : nom de résidu présent dans
  `residues/data/standard_amino_acids.yaml` ;
- **eau** : nom présent dans `residues/data/structure_components.yaml` (`water`) ;
- **ion** : résidu d'un seul atome, qui n'est pas une eau ;
- **résidu de type acide aminé** (pour la numérotation) : acide aminé standard,
  ou résidu portant les atomes de squelette `N`, `CA` et `C` (nomenclature PDB) ;
  c'est le cas de la sélénométhionine `MSE` ou de l'acide aminoisobutyrique `AIB` ;
- **tout le reste** est un **composant non reconnu**, rapporté par son nom (ligand,
  glycane, nucléotide, résidu modifié, groupe terminal) : il n'est jamais
  classé par supposition.

Un **trou de numérotation** est un saut de plus d'une unité entre deux résidus
consécutifs de type acide aminé d'une même chaîne. Ce n'est ni la preuve d'un
résidu manquant ni celle d'une rupture (numérotation de la chymotrypsine,
protéine de fusion) : c'est le rapport d'extraction
(`pepforge.io.peptide_extraction`) qui mesure la liaison C–N de chaque jonction.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from itertools import pairwise
from typing import Any

from pepforge.core.structure import (
    STRUCTURE_FORMAT,
    STRUCTURE_MODEL_COUNT,
    Atom,
    ResidueId,
    Structure,
)
from pepforge.residues.library import ResidueLibrary, load_default_library
from pepforge.residues.structure_components import (
    StructureComponents,
    load_structure_components,
)

BACKBONE_ATOMS = frozenset({"N", "CA", "C"})
"""Atomes de squelette d'un acide aminé α dans la nomenclature PDB."""

STRUCTURE_SOURCE = "structure_source"
"""Clé de `Peptide.metadata` : fichier dont le peptide a été tiré."""

STRUCTURE_CHAIN = "structure_chain"
"""Clé de `Peptide.metadata` : chaîne dont le peptide a été tiré."""

_FROM_STRUCTURE_DEFAULT = "absente de la structure : défaut de la bibliothèque"


@dataclass(frozen=True)
class AltLocSite:
    """Résidu qui porte des conformations alternatives."""

    residue: ResidueId
    alt_ids: tuple[str, ...]
    atoms: int
    """Nombre d'atomes du résidu marqués d'une conformation alternative."""

    residue_names: tuple[str, ...] = ()
    """Noms de résidu différents selon la conformation (microhétérogénéité) ; vide sinon."""


@dataclass(frozen=True)
class ChainInventory:
    """Contenu d'une chaîne."""

    chain_id: str
    atoms: int
    standard_residues: int
    amino_acid_like_residues: int
    waters: int
    ions: tuple[str, ...]
    """Ions, sous la forme `A:ZN701`."""

    unrecognized: dict[str, int]
    """Composants non reconnus de la chaîne, par nom, avec leur nombre."""

    first_residue: str | None
    last_residue: str | None
    numbering_gaps: tuple[tuple[str, str], ...]
    """Couples de résidus consécutifs séparés par un saut de numérotation."""

    insertion_codes: tuple[str, ...]


@dataclass(frozen=True)
class StructureInventory:
    """Inventaire d'une structure, et décisions qu'elle demande."""

    source: str
    format: str | None
    model: int | None
    model_count: int | None
    atoms: int
    hydrogens: int
    chains: tuple[ChainInventory, ...]
    alternate_locations: tuple[AltLocSite, ...]
    declared_missing_residues: tuple[str, ...]
    connections: dict[str, int]
    """Liaisons déclarées par le fichier, par nature."""

    atoms_without_occupancy: int
    atoms_with_zero_occupancy: int
    pending_decisions: tuple[str, ...] = field(default=())

    def to_dict(self) -> dict[str, Any]:
        """Forme sérialisable (YAML, rapport)."""
        return asdict(self)


def _is_amino_acid_like(name: str, atoms: list[Atom], library: ResidueLibrary) -> bool:
    return library.has_code(name) or {a.name for a in atoms} >= BACKBONE_ATOMS


def _is_ion(
    rid: ResidueId, atoms: list[Atom], library: ResidueLibrary, components: StructureComponents
) -> bool:
    """Résidu d'un seul atome (toutes conformations confondues), ni acide aminé ni groupe terminal.

    Le groupe amino terminal `NH2` n'a qu'un atome lourd : sans l'exclusion des
    groupes terminaux, il serait compté comme ion et l'extrémité amidée perdue.
    """
    name = rid.name.upper()
    return (
        len({a.name for a in atoms}) == 1
        and not library.has_code(name)
        and name not in components.n_terminal_caps
        and name not in components.c_terminal_caps
    )


def _residues_by_chain(structure: Structure) -> dict[str, dict[ResidueId, list[Atom]]]:
    chains: dict[str, dict[ResidueId, list[Atom]]] = {}
    for atom in structure.atoms:
        chains.setdefault(atom.chain_id, {}).setdefault(atom.residue_id, []).append(atom)
    return chains


def _alt_sites(structure: Structure) -> tuple[AltLocSite, ...]:
    """Résidus à conformations alternatives, regroupés par position (chaîne, numéro, insertion)."""
    by_position: dict[tuple[str, int, str], list[Atom]] = {}
    for atom in structure.atoms:
        by_position.setdefault(
            (atom.chain_id, atom.residue_number, atom.insertion_code), []
        ).append(atom)
    sites: list[AltLocSite] = []
    for (chain, number, icode), atoms in by_position.items():
        alternates = [a for a in atoms if a.alt_loc]
        if not alternates:
            continue
        names = tuple(dict.fromkeys(a.residue_name for a in atoms))
        sites.append(
            AltLocSite(
                residue=ResidueId(chain, number, icode, "/".join(names)),
                alt_ids=tuple(sorted({a.alt_loc for a in alternates if a.alt_loc})),
                atoms=len(alternates),
                residue_names=names if len(names) > 1 else (),
            )
        )
    return tuple(sites)


def _chain_inventory(
    chain_id: str,
    residues: dict[ResidueId, list[Atom]],
    library: ResidueLibrary,
    components: StructureComponents,
) -> ChainInventory:
    waters = 0
    ions: list[str] = []
    unrecognized: Counter[str] = Counter()
    standard = 0
    amino_like: list[ResidueId] = []
    seen_positions: set[tuple[int, str]] = set()
    for rid, atoms in residues.items():
        position = (rid.number, rid.insertion_code)
        if rid.name.upper() in components.water:
            waters += 1
            continue
        if _is_ion(rid, atoms, library, components):
            ions.append(rid.label())
            continue
        if library.has_code(rid.name):
            standard += 1
        else:
            unrecognized[rid.name] += 1
        if _is_amino_acid_like(rid.name, atoms, library) and position not in seen_positions:
            seen_positions.add(position)
            amino_like.append(rid)
    gaps = tuple(
        (before.label(), after.label())
        for before, after in pairwise(amino_like)
        if after.number - before.number > 1 or (after.number < before.number)
    )
    insertions = tuple(rid.label() for rid in residues if rid.insertion_code)
    return ChainInventory(
        chain_id=chain_id,
        atoms=sum(len(a) for a in residues.values()),
        standard_residues=standard,
        amino_acid_like_residues=len(amino_like),
        waters=waters,
        ions=tuple(ions),
        unrecognized=dict(sorted(unrecognized.items())),
        first_residue=amino_like[0].label() if amino_like else None,
        last_residue=amino_like[-1].label() if amino_like else None,
        numbering_gaps=gaps,
        insertion_codes=insertions,
    )


def _decisions(inv: StructureInventory) -> tuple[str, ...]:
    decisions: list[str] = []
    if inv.model_count and inv.model_count > 1:
        decisions.append(
            f"Fichier à {inv.model_count} modèles : seul le modèle {inv.model} est lu ; "
            "vérifier que c'est le modèle voulu."
        )
    if inv.alternate_locations:
        micro = [s.residue.label() for s in inv.alternate_locations if s.residue_names]
        text = (
            f"{len(inv.alternate_locations)} résidu(s) à conformations alternatives : "
            "en choisir une par résidu avant toute paramétrisation."
        )
        if micro:
            text += f" Résidus dont le nom change selon la conformation : {', '.join(micro)}."
        decisions.append(text)
    if inv.declared_missing_residues:
        decisions.append(
            f"{len(inv.declared_missing_residues)} résidu(s) déclaré(s) non observé(s) par le "
            "fichier : les reconstruire, ou tronquer et coiffer les extrémités."
        )
    gaps = [g for c in inv.chains for g in c.numbering_gaps]
    if gaps:
        decisions.append(
            f"{len(gaps)} saut(s) de numérotation, qui ne prouvent pas une rupture : le "
            "rapport d'extraction (--chain) mesure la liaison C–N de chaque jonction "
            f"({'; '.join(f'{a} → {b}' for a, b in gaps[:5])}{' ...' if len(gaps) > 5 else ''})."
        )
    unrecognized = Counter[str]()
    for chain in inv.chains:
        unrecognized.update(chain.unrecognized)
    if unrecognized:
        names = ", ".join(f"{k}×{v}" for k, v in sorted(unrecognized.items()))
        decisions.append(
            f"Composants non reconnus par la bibliothèque : {names} ; pour chacun, "
            "le garder (et le paramétrer) ou le retirer."
        )
    ions = [ion for c in inv.chains for ion in c.ions]
    if ions:
        decisions.append(f"{len(ions)} ion(s) : les garder ou les retirer ({', '.join(ions[:8])}).")
    waters = sum(c.waters for c in inv.chains)
    if waters:
        decisions.append(
            f"{waters} molécule(s) d'eau cristallographique : les garder ou les retirer."
        )
    if inv.hydrogens == 0:
        decisions.append(
            "Aucun hydrogène : ils devront être ajoutés selon les états de protonation."
        )
    else:
        decisions.append(
            f"{inv.hydrogens} hydrogène(s) présents : vérifier qu'ils correspondent aux états "
            "de protonation choisis, ou les retirer pour les reconstruire."
        )
    if inv.atoms_with_zero_occupancy:
        decisions.append(
            f"{inv.atoms_with_zero_occupancy} atome(s) d'occupation nulle "
            "(position non déterminée)."
        )
    if inv.atoms_without_occupancy:
        decisions.append(f"{inv.atoms_without_occupancy} atome(s) sans occupation dans le fichier.")
    if inv.connections:
        detail = ", ".join(f"{k} ×{v}" for k, v in sorted(inv.connections.items()))
        decisions.append(f"Liaisons déclarées par le fichier ({detail}) : à confirmer.")
    return tuple(decisions)


def inventory(
    structure: Structure,
    library: ResidueLibrary | None = None,
    components: StructureComponents | None = None,
) -> StructureInventory:
    """Inventaire complet de `structure`, sans la modifier."""
    library = library or load_default_library()
    components = components or load_structure_components(library)
    chains = tuple(
        _chain_inventory(chain_id, residues, library, components)
        for chain_id, residues in _residues_by_chain(structure).items()
    )
    model_count = structure.metadata.get(STRUCTURE_MODEL_COUNT)
    draft = StructureInventory(
        source=structure.source,
        format=structure.metadata.get(STRUCTURE_FORMAT),
        model=structure.model,
        model_count=int(model_count) if model_count else None,
        atoms=len(structure.atoms),
        hydrogens=sum(1 for a in structure.atoms if a.element.upper() in ("H", "D")),
        chains=chains,
        alternate_locations=_alt_sites(structure),
        declared_missing_residues=tuple(r.label() for r in structure.declared_missing_residues),
        connections=dict(sorted(Counter(c.kind for c in structure.connections).items())),
        atoms_without_occupancy=sum(1 for a in structure.atoms if a.occupancy is None),
        atoms_with_zero_occupancy=sum(1 for a in structure.atoms if a.occupancy == 0.0),
    )
    return StructureInventory(**{**draft.__dict__, "pending_decisions": _decisions(draft)})
