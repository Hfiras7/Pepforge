"""Connectivité covalente d'une structure importée : liaisons déclarées et liaisons possibles.

Deux natures d'information, jamais confondues :

- une **liaison déclarée** vient du fichier (LINK, SSBOND, `struct_conn`, CONECT) :
  elle est confrontée aux coordonnées (`assess_declared`) ;
- une **liaison possible** vient de la seule géométrie : deux atomes lourds de
  résidus différents assez proches pour être liés (`find_contacts`). Elle ne
  prouve rien : elle est rapportée pour être soumise à une décision explicite.

Chaque évaluation (`BondAssessment`) porte les deux atomes (résidu d'origine,
chaîne, nom, élément, conformation), la distance mesurée, la source de
l'information et le niveau de certitude qui en découle. Aucun seuil n'est écrit
ici : ils viennent de `residues/data/covalent_geometry.yaml`, avec leurs sources
(`CovalentCriteria`).

Les hydrogènes ne servent jamais à la détection : leur position est souvent
calculée, pas observée. Les eaux non plus. Un métal ne reçoit jamais de liaison
covalente possible : son contact avec un autre atome est rapporté comme contact
métallique, de nature non déterminée.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from enum import StrEnum
from itertools import product

from pepforge.core.structure import Atom, Connection, ConnectionKind, ResidueId
from pepforge.residues.structure_components import CovalentCriteria

HYDROGEN_ELEMENTS = frozenset({"H", "D"})
"""Éléments ignorés par la détection (hydrogène et deutérium)."""

ASSESSED_KINDS = frozenset(
    {ConnectionKind.COVALENT, ConnectionKind.DISULFIDE, ConnectionKind.UNKNOWN}
)
"""Natures de liaison déclarée confrontées aux rayons covalents. Une coordination
métallique ou une liaison hydrogène n'est pas une liaison covalente : elle est
gardée et rapportée, sans mesure contre ces rayons."""


class BondEvidence(StrEnum):
    """Source de l'information sur une liaison."""

    DECLARED = "declared"
    """Déclarée par le fichier (LINK, SSBOND, `struct_conn`, CONECT)."""

    GEOMETRIC = "geometric"
    """Déduite de la seule géométrie : deux atomes assez proches pour être liés."""

    UNDETERMINED = "undetermined"
    """Contact avec un métal : coordination ou liaison, la géométrie ne le dit pas."""


class BondStatus(StrEnum):
    """Verdict de la géométrie sur une liaison déclarée ou possible."""

    ACCEPTED = "accepted"
    """Liaison déclarée, et distance compatible avec une liaison covalente."""

    POSSIBLE = "possible"
    """Aucune déclaration, distance compatible : liaison possible, non prouvée."""

    INCOMPATIBLE = "incompatible"
    """Distance hors de la fenêtre covalente : liaison déclarée trop longue, ou
    recouvrement de deux atomes (plus proches que la distance minimale)."""

    INSUFFICIENT_DATA = "insufficient_data"
    """Mesure impossible : atome absent ou d'occupation nulle, conformations en
    désaccord, liaison vers une copie de symétrie."""

    NOT_ASSESSED = "not_assessed"
    """Coordination métallique ou liaison hydrogène : pas de mesure contre les rayons covalents."""


class Certainty(StrEnum):
    """Niveau de certitude, déduit de la source et du verdict de la géométrie."""

    CONFIRMED = "confirmed"
    """Déclarée par le fichier et compatible avec les coordonnées."""

    DECLARED_UNVERIFIED = "declared_unverified"
    """Déclarée, mais la géométrie ne peut ni la confirmer ni l'infirmer."""

    CONTRADICTED = "contradicted"
    """Déclarée, mais les coordonnées la contredisent (distance incompatible)."""

    POSSIBLE = "possible"
    """Proximité seule : liaison possible, jamais une preuve."""

    OVERLAP = "overlap"
    """Deux atomes qui se recouvrent : ni liaison ni contact normal."""

    UNDETERMINED = "undetermined"
    """Contact métallique : nature non déterminée."""


@dataclass(frozen=True)
class AtomRef:
    """Atome désigné par ses identifiants d'origine, jamais renumérotés."""

    residue: ResidueId
    atom: str
    element: str = ""
    """Élément chimique ; vide quand l'atome est désigné par une déclaration sans être présent."""

    alt_loc: str = ""

    def label(self) -> str:
        """Forme lisible : `A:LYS12:NZ`, ou `A:LYS12:NZ[B]` en conformation B."""
        return f"{self.residue.label()}:{self.atom}" + (f"[{self.alt_loc}]" if self.alt_loc else "")

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.residue.label(),
            "chain": self.residue.chain_id,
            "residue": self.residue.name,
            "number": self.residue.number,
            "insertion_code": self.residue.insertion_code,
            "atom": self.atom,
            "element": self.element,
            "alt_loc": self.alt_loc,
        }


@dataclass(frozen=True)
class BondAssessment:
    """Une liaison déclarée ou possible, avec tout ce qui permet de la vérifier."""

    first: AtomRef
    second: AtomRef
    evidence: BondEvidence
    status: BondStatus
    distance: float | None
    """Distance mesurée (Å) ; la plus courte si plusieurs couples de conformations ont
    été mesurés ; `None` si aucune mesure n'est possible."""

    max_distance: float | None
    """Borne haute de la fenêtre covalente appliquée (Å), `None` si non appliquée."""

    reason: str
    kind: str = ""
    """Nature déclarée (`covalent`, `disulfide`, `unknown`...) ; vide pour une proximité."""

    records: tuple[str, ...] = ()
    """Enregistrements du fichier qui la déclarent ; vide pour une proximité."""

    reported_distance: float | None = None
    """Distance écrite dans le fichier, s'il en donne une."""

    @property
    def certainty(self) -> Certainty:
        if self.evidence is BondEvidence.UNDETERMINED:
            return Certainty.UNDETERMINED
        if self.evidence is BondEvidence.GEOMETRIC:
            return (
                Certainty.OVERLAP if self.status is BondStatus.INCOMPATIBLE else Certainty.POSSIBLE
            )
        if self.status is BondStatus.ACCEPTED:
            return Certainty.CONFIRMED
        if self.status is BondStatus.INCOMPATIBLE:
            return Certainty.CONTRADICTED
        return Certainty.DECLARED_UNVERIFIED

    @property
    def residues(self) -> frozenset[ResidueId]:
        return frozenset({self.first.residue, self.second.residue})

    def label(self) -> str:
        return f"{self.first.label()}–{self.second.label()}"

    def to_dict(self) -> dict[str, object]:
        return {
            "first": self.first.to_dict(),
            "second": self.second.to_dict(),
            "distance": None if self.distance is None else round(self.distance, 3),
            "max_distance": None if self.max_distance is None else round(self.max_distance, 3),
            "evidence": self.evidence.value,
            "records": list(self.records),
            "kind": self.kind,
            "reported_distance": self.reported_distance,
            "status": self.status.value,
            "certainty": self.certainty.value,
            "reason": self.reason,
        }


def involves_metal(
    connection: Connection,
    atoms_by_residue: Mapping[ResidueId, list[Atom]],
    criteria: CovalentCriteria,
) -> bool:
    """Vrai si l'un des deux atomes déclarés est un métal (d'après l'élément du fichier)."""
    for rid, name in (
        (connection.first, connection.first_atom),
        (connection.second, connection.second_atom),
    ):
        if any(
            a.name == name and criteria.is_metal(_element(a)) for a in atoms_by_residue.get(rid, [])
        ):
            return True
    return False


def residue_pair_label(a: ResidueId, b: ResidueId) -> str:
    """Étiquette d'une paire de résidus, indépendante de l'ordre : `A:CYS3–A:CYS10`."""
    first, second = sorted((a, b))
    return f"{first.label()}–{second.label()}"


def _element(atom: Atom) -> str:
    return atom.element.strip().upper()


def is_heavy(atom: Atom) -> bool:
    return _element(atom) not in HYDROGEN_ELEMENTS


def _compatible(a: Atom, b: Atom) -> bool:
    """Deux atomes peuvent coexister : même conformation, ou l'un des deux sans conformation."""
    return not a.alt_loc or not b.alt_loc or a.alt_loc == b.alt_loc


def _distance(a: Atom, b: Atom) -> float:
    return math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))


def _ref(atom: Atom) -> AtomRef:
    return AtomRef(atom.residue_id, atom.name, _element(atom), atom.alt_loc or "")


def _copies(atoms: Iterable[Atom], name: str, alt_loc: str) -> list[Atom]:
    """Copies de l'atome `name` : celle de la conformation déclarée, ou toutes."""
    return [a for a in atoms if a.name == name and (not alt_loc or a.alt_loc == alt_loc)]


def assess_declared(
    connection: Connection,
    atoms_by_residue: Mapping[ResidueId, list[Atom]],
    criteria: CovalentCriteria,
) -> BondAssessment:
    """Confronte une liaison déclarée aux coordonnées du fichier."""
    first_atoms = _copies(
        atoms_by_residue.get(connection.first, []), connection.first_atom, connection.first_alt_loc
    )
    second_atoms = _copies(
        atoms_by_residue.get(connection.second, []),
        connection.second_atom,
        connection.second_alt_loc,
    )

    def ref(rid: ResidueId, name: str, alt: str, found: list[Atom]) -> AtomRef:
        element = _element(found[0]) if found else ""
        return AtomRef(rid, name, element, alt)

    first = ref(connection.first, connection.first_atom, connection.first_alt_loc, first_atoms)
    second = ref(connection.second, connection.second_atom, connection.second_alt_loc, second_atoms)

    def result(
        status: BondStatus, reason: str, distance: float | None, limit: float | None
    ) -> BondAssessment:
        return BondAssessment(
            first,
            second,
            BondEvidence.DECLARED,
            status,
            distance,
            limit,
            reason,
            connection.kind,
            connection.records,
            connection.reported_distance,
        )

    if any(criteria.is_metal(r.element) for r in (first, second)):
        return result(
            BondStatus.NOT_ASSESSED,
            "liaison déclarée avec un métal : coordination ou liaison, non confrontée aux "
            "rayons covalents ; le métal n'est jamais rattaché comme liaison covalente",
            None,
            None,
        )
    if connection.kind not in ASSESSED_KINDS:
        return result(
            BondStatus.NOT_ASSESSED,
            f"liaison déclarée de nature {connection.kind} : non confrontée aux rayons covalents",
            None,
            None,
        )
    if connection.symmetry_mate:
        return result(
            BondStatus.INSUFFICIENT_DATA,
            "liaison vers une copie de symétrie : non mesurable sur les coordonnées du fichier",
            None,
            None,
        )
    missing = [
        r.label() for r, found in ((first, first_atoms), (second, second_atoms)) if not found
    ]
    if missing:
        return result(
            BondStatus.INSUFFICIENT_DATA,
            f"atome(s) déclaré(s) absent(s) des coordonnées : {', '.join(missing)}",
            None,
            None,
        )
    limit = criteria.max_distance(first.element, second.element)
    pairs = [
        (a, b)
        for a, b in product(first_atoms, second_atoms)
        if _compatible(a, b)
        and (a.occupancy is None or a.occupancy > 0)
        and (b.occupancy is None or b.occupancy > 0)
    ]
    if not pairs:
        return result(
            BondStatus.INSUFFICIENT_DATA,
            "aucun couple d'atomes mesurable (conformations incompatibles ou occupation nulle)",
            None,
            limit,
        )
    distances = [_distance(a, b) for a, b in pairs]
    inside = [criteria.min_distance <= d <= limit for d in distances]
    shortest = min(distances)
    window = f"fenêtre {criteria.min_distance}–{limit:.2f} Å"
    if all(inside):
        return result(
            BondStatus.ACCEPTED,
            f"déclarée ({'+'.join(connection.records) or 'programme'}), "
            f"{shortest:.2f} Å dans la {window}",
            shortest,
            limit,
        )
    if not any(inside):
        what = "recouvrement" if max(distances) < criteria.min_distance else "distance trop longue"
        return result(
            BondStatus.INCOMPATIBLE,
            f"déclarée mais {what} : {', '.join(f'{d:.2f}' for d in sorted(distances))} Å "
            f"hors de la {window} (déclaration périmée ou coordonnées erronées)",
            shortest,
            limit,
        )
    return result(
        BondStatus.INSUFFICIENT_DATA,
        f"les conformations se contredisent : {', '.join(f'{d:.2f}' for d in sorted(distances))} Å "
        f"pour la {window}",
        shortest,
        limit,
    )


class _Grid:
    """Recherche de voisins par cellules cubiques, en Python pur, déterministe."""

    def __init__(self, atoms: Iterable[Atom], cell: float) -> None:
        self.cell = cell
        self.cells: dict[tuple[int, int, int], list[Atom]] = {}
        for atom in atoms:
            self.cells.setdefault(self._key(atom), []).append(atom)

    def _key(self, atom: Atom) -> tuple[int, int, int]:
        return (
            math.floor(atom.x / self.cell),
            math.floor(atom.y / self.cell),
            math.floor(atom.z / self.cell),
        )

    def near(self, atom: Atom) -> Iterator[Atom]:
        x, y, z = self._key(atom)
        for dx, dy, dz in product((-1, 0, 1), repeat=3):
            yield from self.cells.get((x + dx, y + dy, z + dz), ())


def find_contacts(
    sources: Iterable[Atom],
    targets: Iterable[Atom],
    criteria: CovalentCriteria,
    *,
    skip: frozenset[tuple[ResidueId, str, ResidueId, str]] = frozenset(),
) -> list[BondAssessment]:
    """Paires d'atomes lourds de résidus différents assez proches pour être liés.

    `sources` sont les atomes à examiner (la chaîne, puis les composants qui s'y
    rattachent), `targets` ceux parmi lesquels chercher. `skip` liste les paires
    attendues (liaison peptidique entre résidus consécutifs, liaison déclarée),
    sous la forme (résidu, atome, résidu, atome), dans les deux sens. Chaque
    paire de résidus et d'atomes n'est rapportée qu'une fois, avec la plus
    courte distance entre conformations compatibles.
    """
    heavy_targets = [a for a in targets if is_heavy(a)]
    grid = _Grid(heavy_targets, criteria.largest_distance)
    best: dict[tuple[ResidueId, str, ResidueId, str], tuple[float, Atom, Atom]] = {}
    for source in sources:
        if not is_heavy(source):
            continue
        for target in grid.near(source):
            if target.residue_id == source.residue_id or not _compatible(source, target):
                continue
            key = (source.residue_id, source.name, target.residue_id, target.name)
            if key in skip:
                continue
            metal = criteria.is_metal(_element(source)) or criteria.is_metal(_element(target))
            limit = (
                criteria.largest_distance
                if metal
                else criteria.max_distance(_element(source), _element(target))
            )
            distance = _distance(source, target)
            if distance > limit:
                continue
            ordered = (
                key
                if (source.residue_id, source.name) <= (target.residue_id, target.name)
                else (target.residue_id, target.name, source.residue_id, source.name)
            )
            if ordered not in best or distance < best[ordered][0]:
                best[ordered] = (distance, source, target)
    found: list[BondAssessment] = []
    for _, (distance, a, b) in sorted(best.items(), key=lambda item: item[0]):
        first, second = sorted((_ref(a), _ref(b)), key=lambda r: (r.residue, r.atom))
        if criteria.is_metal(first.element) or criteria.is_metal(second.element):
            found.append(
                BondAssessment(
                    first,
                    second,
                    BondEvidence.UNDETERMINED,
                    BondStatus.NOT_ASSESSED,
                    distance,
                    None,
                    f"contact métallique à {distance:.2f} Å : coordination ou liaison, "
                    "non déterminé ; le métal n'est jamais rattaché comme liaison covalente",
                )
            )
            continue
        limit = criteria.max_distance(first.element, second.element)
        if distance < criteria.min_distance:
            found.append(
                BondAssessment(
                    first,
                    second,
                    BondEvidence.GEOMETRIC,
                    BondStatus.INCOMPATIBLE,
                    distance,
                    limit,
                    f"recouvrement : {distance:.2f} Å, sous la distance minimale "
                    f"{criteria.min_distance} Å ; ni liaison ni contact normal",
                )
            )
            continue
        found.append(
            BondAssessment(
                first,
                second,
                BondEvidence.GEOMETRIC,
                BondStatus.POSSIBLE,
                distance,
                limit,
                f"aucune liaison déclarée ; {distance:.2f} Å ≤ {limit:.2f} Å "
                f"(somme des rayons covalents + {criteria.tolerance} Å) : liaison possible, "
                "non prouvée",
            )
        )
    return found
