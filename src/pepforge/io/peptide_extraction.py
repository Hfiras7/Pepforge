"""Extraction d'un peptide depuis une chaîne de structure, avec un rapport toujours produit.

Le rapport d'extraction (`ExtractionReport`) classe, sans rien modifier :

- chaque **jonction** entre deux résidus consécutifs de la chaîne, d'après la
  distance C(i)–N(i+1) mesurée dans le fichier et le seuil sourcé de
  `residues/data/structure_geometry.yaml` (`JunctionCategory`) ;
- chaque **composant** qui n'est pas un résidu du squelette : groupe terminal,
  modification covalente d'une chaîne latérale (glycane), ligand non lié, ion,
  eau (`ComponentCategory`) ;
- les **résidus déclarés non observés** (`REMARK 465`, `pdbx_unobs_or_zero_occ_residues`),
  placés aux extrémités ou dans une rupture, et les **atomes absents**, déclarés
  par le fichier ou manquants au squelette, rapportés à part.

Chaque jonction et chaque composant porte la décision appliquée (`Decision`).
Le mode par défaut est strict : toute rupture et toute extrémité non observée
est refusée, sauf décision explicite de l'utilisateur (`ExtractionDecisions`).
Rien n'est jamais reconnecté, reconstruit ni corrigé : `split` découpe la
chaîne en segments continus, `observed_segment` garde le segment observé.

Règles de classement, toutes structurelles ou tirées d'un fichier de données :

- **résidu du squelette** : acide aminé de la bibliothèque, résidu portant les
  atomes `N`, `CA`, `C` (acide aminé modifié comme `MSE`), ou groupe terminal
  de `structure_components.yaml` ;
- **modification de chaîne latérale** : composant relié à un résidu du
  squelette, directement ou par d'autres composants, par une liaison covalente
  **déclarée** par le fichier (`LINK`, `struct_conn`) ;
- **ion** : composant d'un seul atome ; **eau** : nom de `structure_components.yaml` ;
- **ligand non lié** : tout autre composant, sans liaison covalente déclarée.
  Sa plus courte distance à la chaîne est mesurée et rapportée, sans seuil :
  un fichier sans table de liaisons (sortie de docking) ne déclare rien, et
  c'est à l'utilisateur de juger.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import StrEnum
from itertools import pairwise, product
from typing import Any

from pepforge.core.peptide import (
    C_TERMINUS_DEFAULTED,
    N_TERMINUS_DEFAULTED,
    Link,
    LinkKind,
    Peptide,
)
from pepforge.core.structure import Atom, ConnectionKind, ResidueId, Structure
from pepforge.exceptions import ExtractionRefusedError, StructureError
from pepforge.report.journal import Journal
from pepforge.residues.library import ResidueLibrary, load_default_library
from pepforge.residues.structure_components import (
    PeptideBondCriterion,
    StructureComponents,
    load_peptide_bond_criterion,
    load_structure_components,
)

BACKBONE_ATOMS = frozenset({"N", "CA", "C"})
"""Atomes qui font d'un résidu un acide aminé α (nomenclature PDB)."""

REPORTED_BACKBONE_ATOMS = ("N", "CA", "C", "O")
"""Atomes de squelette dont l'absence est rapportée pour chaque résidu du squelette."""

STRUCTURE_SOURCE = "structure_source"
"""Clé de `Peptide.metadata` : fichier dont le peptide a été tiré."""

STRUCTURE_CHAIN = "structure_chain"
"""Clé de `Peptide.metadata` : chaîne dont le peptide a été tiré."""

STRUCTURE_SEGMENT = "structure_segment"
"""Clé de `Peptide.metadata` : premier et dernier résidu du segment extrait (`A:ALA1-A:GLY9`)."""

STRUCTURE_EDGES = "structure_edges"
"""Clé de `Peptide.metadata` : nature des bords N et C du segment (`chemical/break`...)."""

_FROM_STRUCTURE_DEFAULT = "absente de la structure : défaut de la bibliothèque"


class JunctionCategory(StrEnum):
    """Nature d'une jonction entre deux résidus consécutifs de la chaîne."""

    CONTINUOUS = "continuous"
    """Liaison C–N présente, numérotation régulière (n → n+1, sans code d'insertion)."""

    NUMBERING_DISCONTINUITY = "numbering_discontinuity"
    """Liaison C–N présente, numérotation irrégulière (saut, retour en arrière, insertion)."""

    BREAK_DECLARED = "break_declared"
    """Pas de liaison C–N, et des résidus déclarés non observés entre les deux."""

    BREAK_UNDECLARED = "break_undeclared"
    """Pas de liaison C–N, et aucun résidu déclaré non observé entre les deux."""

    CONTINUITY_UNMEASURABLE = "continuity_unmeasurable"
    """Atome `C` ou `N` absent, ou conformations alternatives qui ne concordent pas."""


BREAKS = (JunctionCategory.BREAK_DECLARED, JunctionCategory.BREAK_UNDECLARED)


class ComponentCategory(StrEnum):
    """Nature d'un composant de la chaîne."""

    STANDARD_RESIDUE = "standard_residue"
    MODIFIED_RESIDUE = "modified_residue"
    """Résidu du squelette absent de la bibliothèque (`MSE`, `HYP`)."""

    TERMINAL_CAP = "terminal_cap"
    SIDE_CHAIN_MODIFICATION = "side_chain_modification"
    UNBOUND_LIGAND = "unbound_ligand"
    ION = "ion"
    WATER = "water"


class Decision(StrEnum):
    """Décision appliquée à une jonction, une extrémité ou un composant."""

    ACCEPT = "accept"
    REJECT = "reject"
    SPLIT = "split"
    OBSERVED_SEGMENT = "observed_segment"
    EXCLUDE = "exclude"
    """Composant gardé hors du peptide (ligand non lié, ion, eau), et listé."""


BREAK_DECISIONS = (Decision.REJECT, Decision.SPLIT)
UNOBSERVED_TERMINUS_DECISIONS = (Decision.REJECT, Decision.OBSERVED_SEGMENT)


@dataclass(frozen=True)
class ExtractionDecisions:
    """Décisions de l'utilisateur, par catégorie. Sans décision, tout est refusé."""

    breaks: Decision = Decision.REJECT
    """`reject`, ou `split` : chaque segment continu devient un peptide distinct."""

    unobserved_termini: Decision = Decision.REJECT
    """`reject`, ou `observed_segment` : seul le segment observé est extrait."""

    def __post_init__(self) -> None:
        for name, value, allowed in (
            ("breaks", self.breaks, BREAK_DECISIONS),
            ("unobserved_termini", self.unobserved_termini, UNOBSERVED_TERMINUS_DECISIONS),
        ):
            try:
                decision = Decision(value)
            except ValueError:
                decision = None
            if decision not in allowed:
                raise ValueError(
                    f"Décision {name} = {value!r} : valeurs possibles "
                    f"{', '.join(d.value for d in allowed)}"
                )
            object.__setattr__(self, name, decision)

    def to_dict(self) -> dict[str, str]:
        return {"breaks": self.breaks.value, "unobserved_termini": self.unobserved_termini.value}


@dataclass(frozen=True)
class Measurement:
    """Distance C(i)–N(i+1) mesurée entre deux atomes du fichier."""

    c_atom: str
    """Atome `C` utilisé, sous la forme `A:ARG217:C` (ou `A:ARG217:C[B]` en conformation B)."""

    n_atom: str
    distance: float
    """Distance en Å."""


@dataclass(frozen=True)
class Junction:
    """Jonction entre deux résidus consécutifs du squelette, dans l'ordre du fichier."""

    before: str
    after: str
    measurements: tuple[Measurement, ...]
    """Une mesure par couple de conformations ; vide si un atome manque."""

    numbering_step: int
    """Écart de numérotation (`after - before`) ; 1 pour une numérotation régulière."""

    insertion_code: bool
    """Vrai si l'un des deux résidus porte un code d'insertion."""

    declared_missing_between: tuple[str, ...]
    category: JunctionCategory
    decision: Decision
    evidence: str

    @property
    def distance(self) -> float | None:
        """Distance C–N si elle est unique (une seule mesure), sinon `None`."""
        return self.measurements[0].distance if len(self.measurements) == 1 else None


@dataclass(frozen=True)
class Component:
    """Composant de la chaîne qui n'est pas un résidu standard."""

    residue: str
    name: str
    category: ComponentCategory
    decision: Decision
    evidence: str
    bonds: tuple[str, ...] = ()
    """Liaisons déclarées du composant (`A:ASN53:ND2–A:NAG904:C1 (covalent)`)."""

    min_distance_to_chain: float | None = None
    """Plus courte distance (Å) entre un atome du composant et un atome du squelette."""


@dataclass(frozen=True)
class UnobservedTerminus:
    """Résidus déclarés non observés avant le premier ou après le dernier résidu observé."""

    side: str
    """`N` ou `C`."""

    residues: tuple[str, ...]
    decision: Decision


@dataclass(frozen=True)
class Segment:
    """Segment continu extrait : ses bornes et la nature de chaque bord."""

    first: str
    last: str
    residues: int
    n_edge: str
    """`chemical` (vraie extrémité de la molécule), `unobserved` ou `break`."""

    c_edge: str


@dataclass(frozen=True)
class ExtractionReport:
    """Rapport d'extraction d'une chaîne, produit même quand l'extraction est refusée."""

    source: str
    chain_id: str
    model: int | None
    criterion: dict[str, Any]
    decisions: dict[str, str]
    residues: tuple[str, ...]
    """Résidus du squelette, groupes terminaux compris, dans l'ordre du fichier."""

    junctions: tuple[Junction, ...]
    components: tuple[Component, ...]
    unobserved_termini: tuple[UnobservedTerminus, ...]
    declared_missing_residues: tuple[str, ...]
    declared_missing_atoms: tuple[str, ...]
    absent_backbone_atoms: tuple[str, ...]
    """Atomes `N`, `CA`, `C`, `O` absents d'un résidu du squelette, déclarés ou non."""

    segments: tuple[Segment, ...] = ()
    refusals: tuple[str, ...] = field(default=())

    @property
    def accepted(self) -> bool:
        return not self.refusals

    def junctions_by_category(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for junction in self.junctions:
            counts[junction.category.value] = counts.get(junction.category.value, 0) + 1
        return dict(sorted(counts.items()))

    def to_dict(self) -> dict[str, Any]:
        """Forme sérialisable (YAML, rapport). Les distances sont arrondies au millième d'Å."""
        data = asdict(self)
        for junction, raw in zip(self.junctions, data["junctions"], strict=True):
            raw["category"] = junction.category.value
            raw["decision"] = junction.decision.value
            raw["distance"] = None if junction.distance is None else round(junction.distance, 3)
            for measure in raw["measurements"]:
                measure["distance"] = round(measure["distance"], 3)
        for component, raw in zip(self.components, data["components"], strict=True):
            raw["category"] = component.category.value
            raw["decision"] = component.decision.value
            if component.min_distance_to_chain is not None:
                raw["min_distance_to_chain"] = round(component.min_distance_to_chain, 3)
        for terminus, raw in zip(self.unobserved_termini, data["unobserved_termini"], strict=True):
            raw["decision"] = terminus.decision.value
        data["junctions_by_category"] = self.junctions_by_category()
        data["accepted"] = self.accepted
        return data


@dataclass(frozen=True)
class ExtractionResult:
    """Peptides extraits (un par segment) et rapport."""

    report: ExtractionReport
    peptides: tuple[Peptide, ...]


# --- Classement -----------------------------------------------------------------------------


def _atom_label(rid: ResidueId, atom: Atom) -> str:
    alt = f"[{atom.alt_loc}]" if atom.alt_loc else ""
    return f"{rid.label()}:{atom.name}{alt}"


def _distance(a: Atom, b: Atom) -> float:
    return math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))


def _is_heavy(atom: Atom) -> bool:
    return atom.element.upper() not in ("H", "D")


def _measure(
    before: tuple[ResidueId, list[Atom]], after: tuple[ResidueId, list[Atom]]
) -> tuple[Measurement, ...]:
    """Toutes les distances C–N entre conformations compatibles (même indicateur, ou sans)."""
    carbons = [a for a in before[1] if a.name == "C"]
    nitrogens = [a for a in after[1] if a.name == "N"]
    return tuple(
        Measurement(_atom_label(before[0], c), _atom_label(after[0], n), _distance(c, n))
        for c, n in product(carbons, nitrogens)
        if not (c.alt_loc and n.alt_loc and c.alt_loc != n.alt_loc)
    )


def _order_key(rid: ResidueId) -> tuple[int, str]:
    return (rid.number, rid.insertion_code)


def _between(missing: ResidueId, before: ResidueId, after: ResidueId) -> bool:
    """Résidu déclaré numéroté strictement entre deux résidus (numéro, puis code d'insertion)."""
    low, high = sorted((_order_key(before), _order_key(after)))
    return low < _order_key(missing) < high


def _junction(
    before: tuple[ResidueId, list[Atom]],
    after: tuple[ResidueId, list[Atom]],
    missing: list[ResidueId],
    criterion: PeptideBondCriterion,
    decisions: ExtractionDecisions,
) -> Junction:
    b, a = before[0], after[0]
    measurements = _measure(before, after)
    step = a.number - b.number
    icode = bool(b.insertion_code or a.insertion_code)
    between = tuple(m.label() for m in missing if _between(m, b, a))
    threshold = criterion.max_c_n_distance
    bonded = [m.distance <= threshold for m in measurements]
    shown = ", ".join(f"{m.distance:.3f} Å" for m in measurements)
    category: JunctionCategory
    if not measurements:
        absent = [
            f"{rid.label()}:{name}"
            for rid, atoms, name in ((b, before[1], "C"), (a, after[1], "N"))
            if not any(x.name == name for x in atoms)
        ]
        category = JunctionCategory.CONTINUITY_UNMEASURABLE
        evidence = f"atome absent : {', '.join(absent)}"
    elif all(bonded):
        regular = step == 1 and not icode
        category = (
            JunctionCategory.CONTINUOUS if regular else JunctionCategory.NUMBERING_DISCONTINUITY
        )
        evidence = f"géométrie : C–N {shown} ≤ {threshold} Å"
        if not regular:
            evidence += f" ; numérotation {step:+d}" + (" avec code d'insertion" if icode else "")
    elif not any(bonded):
        category = JunctionCategory.BREAK_DECLARED if between else JunctionCategory.BREAK_UNDECLARED
        evidence = f"géométrie : C–N {shown} > {threshold} Å ; "
        evidence += (
            f"en-tête : {len(between)} résidu(s) déclaré(s) non observé(s) entre les deux"
            if between
            else "en-tête : aucun résidu déclaré non observé entre les deux"
        )
    else:
        category = JunctionCategory.CONTINUITY_UNMEASURABLE
        evidence = f"conformations alternatives discordantes : C–N {shown}, seuil {threshold} Å"
    if category in (JunctionCategory.CONTINUOUS, JunctionCategory.NUMBERING_DISCONTINUITY):
        decision = Decision.ACCEPT
    elif category is JunctionCategory.CONTINUITY_UNMEASURABLE:
        decision = Decision.REJECT
    else:
        decision = decisions.breaks
    return Junction(
        before=b.label(),
        after=a.label(),
        measurements=measurements,
        numbering_step=step,
        insertion_code=icode,
        declared_missing_between=between,
        category=category,
        decision=decision,
        evidence=evidence,
    )


def _residues_of_chain(structure: Structure, chain_id: str) -> dict[ResidueId, list[Atom]]:
    residues: dict[ResidueId, list[Atom]] = {}
    for atom in structure.atoms:
        if atom.chain_id == chain_id:
            residues.setdefault(atom.residue_id, []).append(atom)
    return residues


def _is_backbone(
    rid: ResidueId, atoms: list[Atom], library: ResidueLibrary, components: StructureComponents
) -> bool:
    name = rid.name.upper()
    return (
        library.has_code(rid.name)
        or name in components.n_terminal_caps
        or name in components.c_terminal_caps
        or {a.name for a in atoms} >= BACKBONE_ATOMS
    )


def _same_position(a: ResidueId, b: ResidueId) -> bool:
    return (a.chain_id, a.number, a.insertion_code) == (b.chain_id, b.number, b.insertion_code)


def _connection_label(conn: Any) -> str:
    return (
        f"{conn.first.label()}:{conn.first_atom}–{conn.second.label()}:{conn.second_atom} "
        f"({conn.kind})"
    )


def _min_distance(atoms: Iterable[Atom], chain_atoms: list[Atom]) -> float | None:
    best: float | None = None
    for atom in atoms:
        for other in chain_atoms:
            d = _distance(atom, other)
            if best is None or d < best:
                best = d
    return best


@dataclass
class _Analysis:
    """État intermédiaire de l'analyse d'une chaîne."""

    backbone: list[tuple[ResidueId, list[Atom]]]
    report: ExtractionReport


def _analyze(
    structure: Structure,
    chain_id: str,
    decisions: ExtractionDecisions,
    library: ResidueLibrary,
    components: StructureComponents,
    criterion: PeptideBondCriterion,
) -> _Analysis:
    residues = _residues_of_chain(structure, chain_id)
    if not residues:
        raise StructureError(
            f"Chaîne {chain_id!r} absente ; chaînes présentes : {structure.chain_ids}"
        )
    refusals: list[str] = []

    # Microhétérogénéité : deux résidus différents à la même position.
    by_position: dict[tuple[int, str], list[str]] = {}
    for rid in residues:
        by_position.setdefault((rid.number, rid.insertion_code), []).append(rid.name)
    for (number, icode), names in by_position.items():
        if len(names) > 1:
            refusals.append(
                f"{chain_id}:{number}{icode} porte plusieurs résidus selon la conformation "
                f"({'/'.join(names)}) : microhétérogénéité, à résoudre avant l'extraction"
            )

    backbone: list[tuple[ResidueId, list[Atom]]] = []
    others: list[tuple[ResidueId, list[Atom]]] = []
    for rid, atoms in residues.items():
        name = rid.name.upper()
        is_water = name in components.water
        if not is_water and _is_backbone(rid, atoms, library, components):
            if not any(_same_position(rid, kept) for kept, _ in backbone):
                backbone.append((rid, atoms))
        else:
            others.append((rid, atoms))
    if not backbone:
        names = sorted({rid.name for rid, _ in others})
        raise StructureError(
            f"Chaîne {chain_id} : aucun acide aminé (ni N, CA, C, ni groupe terminal) ; "
            f"chaîne non peptidique, composants : {', '.join(names[:8])}"
        )
    backbone_ids = {rid for rid, _ in backbone}

    # Groupes terminaux : seulement à leur extrémité.
    for position, (rid, _) in enumerate(backbone):
        name = rid.name.upper()
        if name in components.n_terminal_caps and position != 0:
            refusals.append(f"groupe terminal N {rid.label()} ailleurs qu'en tête de chaîne")
        if name in components.c_terminal_caps and position != len(backbone) - 1:
            refusals.append(f"groupe terminal C {rid.label()} ailleurs qu'en fin de chaîne")

    # Jonctions.
    missing = [m for m in structure.declared_missing_residues if m.chain_id == chain_id]
    junctions = tuple(
        _junction(before, after, missing, criterion, decisions)
        for before, after in pairwise(backbone)
    )

    # Résidus déclarés non observés : extrémités, ruptures, ou non placés.
    first, last = backbone[0][0], backbone[-1][0]
    n_side = [m for m in missing if _order_key(m) < _order_key(first)]
    c_side = [m for m in missing if _order_key(m) > _order_key(last)]
    unobserved: list[UnobservedTerminus] = []
    for side, group in (("N", n_side), ("C", c_side)):
        if group:
            unobserved.append(
                UnobservedTerminus(
                    side, tuple(m.label() for m in group), decisions.unobserved_termini
                )
            )
    placed_in_breaks = {
        label
        for j in junctions
        if j.category is JunctionCategory.BREAK_DECLARED
        for label in j.declared_missing_between
    }
    placed = placed_in_breaks | {m.label() for m in n_side} | {m.label() for m in c_side}
    for m in missing:
        if m.label() not in placed:
            refusals.append(
                f"résidu déclaré non observé {m.label()} sans place dans la chaîne observée "
                "(ni à une extrémité, ni dans une rupture mesurée) : en-tête et coordonnées "
                "se contredisent"
            )

    # Composants.
    covalent = [
        c
        for c in structure.connections
        if c.kind in (ConnectionKind.COVALENT, ConnectionKind.DISULFIDE, ConnectionKind.UNKNOWN)
    ]
    other_ids = {rid for rid, _ in others}
    attached: set[ResidueId] = set()
    frontier = set(backbone_ids)
    while frontier:
        found: set[ResidueId] = set()
        for conn in covalent:
            for here, there in ((conn.first, conn.second), (conn.second, conn.first)):
                if here in frontier and there in other_ids and there not in attached:
                    found.add(there)
        attached |= found
        frontier = found
    backbone_atoms = [a for _, atoms in backbone for a in atoms if _is_heavy(a)]
    component_list: list[Component] = []
    for rid, atoms in others:
        name = rid.name.upper()
        bonds = tuple(
            _connection_label(c) for c in structure.connections if rid in (c.first, c.second)
        )
        if name in components.water:
            component_list.append(
                Component(
                    rid.label(),
                    rid.name,
                    ComponentCategory.WATER,
                    Decision.EXCLUDE,
                    "nom d'eau (structure_components.yaml)",
                    bonds,
                )
            )
        elif rid in attached:
            component_list.append(
                Component(
                    rid.label(),
                    rid.name,
                    ComponentCategory.SIDE_CHAIN_MODIFICATION,
                    Decision.REJECT,
                    "lié à la chaîne par une liaison covalente déclarée : modification "
                    "non représentable dans le peptide",
                    bonds,
                )
            )
            refusals.append(
                f"modification covalente {rid.label()} liée à la chaîne ({'; '.join(bonds)}) : "
                "non représentable dans le peptide, jamais retirée en silence"
            )
        else:
            single = len({a.name for a in atoms}) == 1
            category = ComponentCategory.ION if single else ComponentCategory.UNBOUND_LIGAND
            evidence = (
                "un seul atome" if single else "plusieurs atomes"
            ) + " ; aucune liaison covalente déclarée avec la chaîne"
            component_list.append(
                Component(
                    rid.label(),
                    rid.name,
                    category,
                    Decision.EXCLUDE,
                    evidence,
                    bonds,
                    _min_distance((a for a in atoms if _is_heavy(a)), backbone_atoms),
                )
            )
    for rid, _ in backbone:
        name = rid.name.upper()
        if name in components.n_terminal_caps or name in components.c_terminal_caps:
            component_list.append(
                Component(
                    rid.label(),
                    rid.name,
                    ComponentCategory.TERMINAL_CAP,
                    Decision.ACCEPT,
                    "groupe terminal (structure_components.yaml)",
                )
            )
        elif not library.has_code(rid.name):
            component_list.append(
                Component(
                    rid.label(),
                    rid.name,
                    ComponentCategory.MODIFIED_RESIDUE,
                    Decision.REJECT,
                    "atomes N, CA, C présents ; résidu absent de la bibliothèque",
                )
            )
            refusals.append(
                f"résidu modifié {rid.label()} dans la chaîne : absent de la bibliothèque"
            )

    # Liaisons déclarées entre résidus du squelette : peptidiques, disulfures, ou refus.
    order = {rid: i for i, (rid, _) in enumerate(backbone)}
    for conn in covalent:
        ends_here = [r for r in (conn.first, conn.second) if r.chain_id == chain_id]
        if not ends_here:
            continue
        i, j = order.get(conn.first), order.get(conn.second)
        if i is None and j is None:
            continue  # entre composants : déjà rapportée avec eux
        if (i is None) != (j is None):
            partner = conn.second if i is not None else conn.first
            if partner in attached:
                continue  # modification de chaîne latérale, rapportée avec le composant
            refusals.append(
                f"liaison déclarée {_connection_label(conn)} vers un résidu hors de la chaîne : "
                "non représentable dans le peptide"
            )
            continue
        assert i is not None and j is not None
        if conn.kind == ConnectionKind.DISULFIDE:
            continue
        (lo, lo_atom), (hi, hi_atom) = sorted(((i, conn.first_atom), (j, conn.second_atom)))
        if hi - lo == 1 and (lo_atom, hi_atom) == ("C", "N"):
            continue  # liaison peptidique déclarée (LINK autour d'un résidu HETATM)
        refusals.append(
            f"liaison déclarée {_connection_label(conn)} hors du squelette : "
            "non représentable dans le peptide"
        )

    # Refus par jonction et par extrémité.
    for junction in junctions:
        if junction.decision is Decision.REJECT:
            refusals.append(
                f"jonction {junction.before} → {junction.after} : {junction.category.value} "
                f"({junction.evidence})"
            )
    for terminus in unobserved:
        if terminus.decision is Decision.REJECT:
            refusals.append(
                f"{len(terminus.residues)} résidu(s) non observé(s) à l'extrémité "
                f"{terminus.side} ({', '.join(terminus.residues[:5])}"
                f"{' ...' if len(terminus.residues) > 5 else ''}) : l'extrémité observée "
                "n'est pas celle de la molécule"
            )

    declared_atoms = tuple(
        m.label() for m in structure.declared_missing_atoms if m.residue.chain_id == chain_id
    )
    absent = tuple(
        f"{rid.label()}:{name}"
        for rid, atoms in backbone
        if not (
            rid.name.upper() in components.n_terminal_caps
            or rid.name.upper() in components.c_terminal_caps
        )
        for name in REPORTED_BACKBONE_ATOMS
        if not any(a.name == name for a in atoms)
    )

    # Segments : découpage aux ruptures acceptées en `split`.
    segments: list[Segment] = []
    start = 0
    n_edge = "unobserved" if n_side else "chemical"
    for k, junction in enumerate(junctions):
        if junction.category in BREAKS:
            segments.append(
                Segment(
                    backbone[start][0].label(),
                    backbone[k][0].label(),
                    k + 1 - start,
                    n_edge,
                    "break",
                )
            )
            start, n_edge = k + 1, "break"
    segments.append(
        Segment(
            backbone[start][0].label(),
            backbone[-1][0].label(),
            len(backbone) - start,
            n_edge,
            "unobserved" if c_side else "chemical",
        )
    )

    report = ExtractionReport(
        source=structure.source,
        chain_id=chain_id,
        model=structure.model,
        criterion={
            "max_c_n_distance": criterion.max_c_n_distance,
            "source": criterion.source,
            "reference_length": criterion.reference_length,
            "reference_source": criterion.reference_source,
        },
        decisions=decisions.to_dict(),
        residues=tuple(rid.label() for rid, _ in backbone),
        junctions=junctions,
        components=tuple(component_list),
        unobserved_termini=tuple(unobserved),
        declared_missing_residues=tuple(m.label() for m in missing),
        declared_missing_atoms=declared_atoms,
        absent_backbone_atoms=absent,
        segments=tuple(segments),
        refusals=tuple(dict.fromkeys(refusals)),
    )
    return _Analysis(backbone, report)


def analyze_chain(
    structure: Structure,
    chain_id: str,
    *,
    decisions: ExtractionDecisions | None = None,
    library: ResidueLibrary | None = None,
    components: StructureComponents | None = None,
    criterion: PeptideBondCriterion | None = None,
) -> ExtractionReport:
    """Rapport d'extraction de la chaîne, sans extraire ni refuser (sauf chaîne absente)."""
    library = library or load_default_library()
    return _analyze(
        structure,
        chain_id,
        decisions or ExtractionDecisions(),
        library,
        components or load_structure_components(library),
        criterion or load_peptide_bond_criterion(),
    ).report


# --- Extraction -----------------------------------------------------------------------------


def _segment_bounds(analysis: _Analysis) -> list[tuple[int, int]]:
    """Indices (début, fin exclue) de chaque segment dans le squelette."""
    bounds: list[tuple[int, int]] = []
    start = 0
    for k, junction in enumerate(analysis.report.junctions):
        if junction.category in BREAKS:
            bounds.append((start, k + 1))
            start = k + 1
    bounds.append((start, len(analysis.backbone)))
    return bounds


def _peptide(
    structure: Structure,
    analysis: _Analysis,
    lo: int,
    hi: int,
    segment: Segment,
    library: ResidueLibrary,
    components: StructureComponents,
) -> Peptide:
    chain_id = analysis.report.chain_id
    where = f"{structure.source or 'structure'}, chaîne {chain_id}"
    part = analysis.backbone[lo:hi]
    n_key: str | None = None
    c_key: str | None = None
    if part[0][0].name.upper() in components.n_terminal_caps:
        n_key = components.n_terminal_caps[part.pop(0)[0].name.upper()]
    if part and part[-1][0].name.upper() in components.c_terminal_caps:
        c_key = components.c_terminal_caps[part.pop()[0].name.upper()]
    if not part:
        raise StructureError(f"{where} : aucun résidu entre les groupes terminaux")
    # Bord d'un segment (rupture, résidus non observés au-delà) : extrémité non
    # chimique, jamais une amine ou un acide supposés, même si un OXT ou un H1 y figure.
    n_open = n_key is None and segment.n_edge != "chemical"
    c_open = c_key is None and segment.c_edge != "chemical"
    # Ailleurs, OXT prouve un acide libre et H1/H2/H3 une amine libre.
    if c_key is None and {a.name for a in part[-1][1]} & components.free_c_terminus_atoms:
        c_key = components.free_c_terminus
    if n_key is None and {a.name for a in part[0][1]} & components.free_n_terminus_hydrogens:
        n_key = components.free_n_terminus
    index = {rid: position for position, (rid, _) in enumerate(part, start=1)}
    links = [
        Link(LinkKind.DISULFIDE, index[conn.first], index[conn.second])
        for conn in structure.connections
        if conn.kind == ConnectionKind.DISULFIDE and conn.first in index and conn.second in index
    ]
    metadata = {
        STRUCTURE_SOURCE: structure.source,
        STRUCTURE_CHAIN: chain_id,
        STRUCTURE_SEGMENT: f"{segment.first}-{segment.last}",
        STRUCTURE_EDGES: f"{segment.n_edge}/{segment.c_edge}",
    }
    if n_key is None and not n_open:
        metadata[N_TERMINUS_DEFAULTED] = _FROM_STRUCTURE_DEFAULT
    if c_key is None and not c_open:
        metadata[C_TERMINUS_DEFAULTED] = _FROM_STRUCTURE_DEFAULT
    n_terminus = (
        library.unobserved_terminus
        if n_open
        else library.n_termini[n_key]
        if n_key
        else library.n_terminus()
    )
    c_terminus = (
        library.unobserved_terminus
        if c_open
        else library.c_termini[c_key]
        if c_key
        else library.c_terminus()
    )
    return Peptide(
        residues=tuple(library.by_code(rid.name) for rid, _ in part),
        n_terminus=n_terminus,
        c_terminus=c_terminus,
        links=tuple(links),
        metadata=metadata,
    )


def _crossing_disulfides(structure: Structure, analysis: _Analysis) -> list[str]:
    """Ponts disulfures de la chaîne dont un bout sort de son segment (autre segment ou chaîne)."""
    segment_of: dict[ResidueId, int] = {}
    for number, (lo, hi) in enumerate(_segment_bounds(analysis)):
        for rid, _ in analysis.backbone[lo:hi]:
            segment_of[rid] = number
    crossing: list[str] = []
    for conn in structure.connections:
        if conn.kind != ConnectionKind.DISULFIDE:
            continue
        ends = (segment_of.get(conn.first), segment_of.get(conn.second))
        if ends == (None, None) or ends[0] == ends[1]:
            continue
        if analysis.report.chain_id in (conn.first.chain_id, conn.second.chain_id):
            crossing.append(
                f"pont disulfure {conn.first.label()}–{conn.second.label()} entre deux "
                "segments ou vers une autre chaîne : non représentable dans un peptide"
            )
    return crossing


def extract_peptides(
    structure: Structure,
    chain_id: str,
    *,
    journal: Journal,
    decisions: ExtractionDecisions | None = None,
    library: ResidueLibrary | None = None,
    components: StructureComponents | None = None,
    criterion: PeptideBondCriterion | None = None,
) -> ExtractionResult:
    """Peptide(s) de la chaîne `chain_id` selon `decisions`, ou `ExtractionRefusedError`.

    Le rapport est inscrit au journal dans les deux cas, et porté par
    l'exception en cas de refus (`ExtractionRefusedError.report`).
    """
    decisions = decisions or ExtractionDecisions()
    library = library or load_default_library()
    components = components or load_structure_components(library)
    criterion = criterion or load_peptide_bond_criterion()
    analysis = _analyze(structure, chain_id, decisions, library, components, criterion)
    crossing = _crossing_disulfides(structure, analysis)
    if crossing:
        report = analysis.report
        analysis.report = ExtractionReport(
            **{**report.__dict__, "refusals": tuple(dict.fromkeys((*report.refusals, *crossing)))}
        )
    report = analysis.report
    if not report.accepted:
        journal.record(
            "extract_peptide_refused",
            "extraction refusée : le peptide lu ne serait pas celui de la molécule",
            **_journal_fields(report),
        )
        where = f"{structure.source or 'structure'}, chaîne {chain_id}"
        raise ExtractionRefusedError(
            f"{where} : extraction refusée ({len(report.refusals)} motif(s)) : "
            + " ; ".join(report.refusals[:3])
            + (" ..." if len(report.refusals) > 3 else ""),
            report,
        )
    peptides = tuple(
        _peptide(structure, analysis, lo, hi, segment, library, components)
        for (lo, hi), segment in zip(_segment_bounds(analysis), report.segments, strict=True)
    )
    journal.record(
        "extract_peptide",
        "correspondance entre une chaîne de la structure et un ou plusieurs peptides",
        **_journal_fields(report),
        peptides=[
            {
                "segment": p.metadata[STRUCTURE_SEGMENT],
                "length": len(p),
                "n_terminus": p.n_terminus.key,
                "edges": p.metadata[STRUCTURE_EDGES],
                "n_terminus_from_structure": N_TERMINUS_DEFAULTED not in p.metadata
                and p.n_terminus.chemical,
                "c_terminus": p.c_terminus.key,
                "c_terminus_from_structure": C_TERMINUS_DEFAULTED not in p.metadata
                and p.c_terminus.chemical,
                "disulfides": [link.positions for link in p.links],
            }
            for p in peptides
        ],
    )
    return ExtractionResult(report, peptides)


def _journal_fields(report: ExtractionReport) -> dict[str, Any]:
    data = report.to_dict()
    return {
        "source": report.source,
        "chain": report.chain_id,
        "model": report.model,
        "criterion": data["criterion"],
        "decisions": data["decisions"],
        "junctions": data["junctions"],
        "components": [c for c in data["components"] if c["category"] != "water"],
        "waters": sum(1 for c in report.components if c.category is ComponentCategory.WATER),
        "unobserved_termini": data["unobserved_termini"],
        "declared_missing_residues": len(report.declared_missing_residues),
        "declared_missing_atoms": list(report.declared_missing_atoms),
        "absent_backbone_atoms": list(report.absent_backbone_atoms),
        "segments": data["segments"],
        "refusals": list(report.refusals),
    }


def structure_to_peptide(
    structure: Structure,
    chain_id: str,
    *,
    journal: Journal,
    library: ResidueLibrary | None = None,
    components: StructureComponents | None = None,
    criterion: PeptideBondCriterion | None = None,
) -> Peptide:
    """Peptide unique de la chaîne, en mode strict (aucune décision) : `extract_peptides`."""
    result = extract_peptides(
        structure,
        chain_id,
        journal=journal,
        library=library,
        components=components,
        criterion=criterion,
    )
    (peptide,) = result.peptides
    return peptide
