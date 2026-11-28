"""Composants chimiques reconnus à la lecture d'une structure (eau, groupes terminaux).

Les données viennent de `data/structure_components.yaml`, avec leur source
(dictionnaire des composants chimiques de la wwPDB). Le fichier est vérifié au
chargement : chaque groupe terminal renvoie à une extrémité connue de
`termini.yaml`, sinon `ReferenceDataError`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pepforge.exceptions import ReferenceDataError
from pepforge.residues.library import ResidueLibrary, load_default_library, read_data_file

DATA_FILE = "structure_components.yaml"


@dataclass(frozen=True)
class StructureComponents:
    """Noms de composants que l'inventaire et la correspondance savent classer."""

    water: frozenset[str]
    n_terminal_caps: dict[str, str]
    """Composant (`ACE`) → clé d'extrémité N (`acetyl`)."""

    c_terminal_caps: dict[str, str]
    """Composant (`NH2`, `NME`) → clé d'extrémité C (`amide`, `n_methylamide`)."""

    free_c_terminus: str
    """Extrémité C prouvée par un atome de `free_c_terminus_atoms` (carboxyle libre)."""

    free_c_terminus_atoms: frozenset[str]
    free_n_terminus: str
    """Extrémité N prouvée par un hydrogène de `free_n_terminus_hydrogens` (amine libre)."""

    free_n_terminus_hydrogens: frozenset[str]
    source: str


def _names(data: dict[str, Any], key: str) -> frozenset[str]:
    value = data.get(key)
    if not isinstance(value, list) or not value:
        raise ReferenceDataError(f"{DATA_FILE} : liste « {key} » absente ou vide")
    return frozenset(str(v).upper() for v in value)


def _caps(data: dict[str, Any], key: str, known: dict[str, Any]) -> dict[str, str]:
    caps: dict[str, str] = {}
    for entry in data.get(key) or []:
        component = str(entry["component"]).upper()
        terminus = str(entry["terminus"])
        if terminus not in known:
            raise ReferenceDataError(
                f"{DATA_FILE} : le groupe {component} renvoie à l'extrémité inconnue {terminus!r}"
            )
        caps[component] = terminus
    return caps


def _terminus(data: dict[str, Any], key: str, known: dict[str, Any]) -> str:
    value = data.get(key)
    if value not in known:
        raise ReferenceDataError(
            f"{DATA_FILE} : « {key} » renvoie à l'extrémité inconnue {value!r}"
        )
    return str(value)


def load_structure_components(library: ResidueLibrary | None = None) -> StructureComponents:
    """Charge et vérifie `structure_components.yaml`."""
    library = library or load_default_library()
    data = read_data_file(DATA_FILE)
    source = data.get("source")
    if not source:
        raise ReferenceDataError(f"{DATA_FILE} : source non déclarée")
    return StructureComponents(
        water=_names(data, "water"),
        n_terminal_caps=_caps(data, "n_terminal_caps", dict(library.n_termini)),
        c_terminal_caps=_caps(data, "c_terminal_caps", dict(library.c_termini)),
        free_c_terminus=_terminus(data, "free_c_terminus", dict(library.c_termini)),
        free_c_terminus_atoms=_names(data, "free_c_terminus_atoms"),
        free_n_terminus=_terminus(data, "free_n_terminus", dict(library.n_termini)),
        free_n_terminus_hydrogens=_names(data, "free_n_terminus_hydrogens"),
        source=str(source),
    )


GEOMETRY_FILE = "structure_geometry.yaml"


@dataclass(frozen=True)
class PeptideBondCriterion:
    """Seuil C(i)–N(i+1) au-delà duquel deux résidus consécutifs ne sont pas liés."""

    max_c_n_distance: float
    """Seuil en Å."""

    source: str
    reference_length: float
    reference_source: str

    def describe(self) -> str:
        return f"C–N ≤ {self.max_c_n_distance} Å"


def load_peptide_bond_criterion(max_c_n_distance: float | None = None) -> PeptideBondCriterion:
    """Charge le seuil de `structure_geometry.yaml` ; une autre valeur peut être imposée.

    Une valeur imposée remplace celle du fichier et sa source devient « valeur
    imposée par l'utilisateur » : elle est recopiée telle quelle dans le rapport.
    """
    data = read_data_file(GEOMETRY_FILE).get("peptide_bond")
    if not isinstance(data, dict):
        raise ReferenceDataError(f"{GEOMETRY_FILE} : section « peptide_bond » absente")
    for key in ("max_c_n_distance", "source", "reference_length", "reference_source"):
        if data.get(key) in (None, ""):
            raise ReferenceDataError(f"{GEOMETRY_FILE} : « peptide_bond.{key} » absent")
    threshold = float(data["max_c_n_distance"])
    source = str(data["source"])
    if max_c_n_distance is not None:
        threshold = float(max_c_n_distance)
        source = f"valeur imposée par l'utilisateur (fichier : {data['max_c_n_distance']} Å)"
    reference = float(data["reference_length"])
    if not threshold > reference:
        raise ReferenceDataError(
            f"Seuil C–N de {threshold} Å inférieur ou égal à la longueur de référence "
            f"{reference} Å : aucune liaison réelle ne serait reconnue"
        )
    return PeptideBondCriterion(threshold, source, reference, str(data["reference_source"]))


COVALENT_FILE = "covalent_geometry.yaml"


@dataclass(frozen=True)
class CovalentCriteria:
    """Fenêtre de distance d'une liaison covalente entre deux atomes lourds (étape C).

    Une liaison est possible entre deux atomes quand
    `min_distance ≤ d ≤ r1 + r2 + tolerance`. Ce critère ne prouve jamais une
    liaison : il dit seulement quand poser la question.
    """

    radii: dict[str, float]
    """Rayon covalent par élément (symbole en majuscules), en Å."""

    tolerance: float
    min_distance: float
    metals: frozenset[str]
    radii_source: str
    tolerance_source: str
    min_distance_source: str
    metals_source: str

    def is_metal(self, element: str) -> bool:
        return element.upper() in self.metals

    def max_distance(self, first: str, second: str) -> float:
        """Distance au-delà de laquelle deux atomes ne sont pas liés, en Å.

        Un élément sans rayon sourcé lève `ReferenceDataError` : jamais de valeur supposée.
        """
        missing = [e for e in (first, second) if e.upper() not in self.radii]
        if missing:
            raise ReferenceDataError(
                f"{COVALENT_FILE} : aucun rayon covalent sourcé pour l'élément "
                f"{', '.join(sorted(set(missing)))} ; la connectivité ne peut pas être vérifiée"
            )
        return self.radii[first.upper()] + self.radii[second.upper()] + self.tolerance

    @property
    def largest_distance(self) -> float:
        """Plus grande distance de liaison possible entre deux éléments connus."""
        return 2 * max(self.radii.values()) + self.tolerance

    def describe(self) -> dict[str, object]:
        """Valeurs et sources, telles qu'elles sont recopiées dans le rapport."""
        return {
            "rule": "min_distance <= d <= r1 + r2 + tolerance",
            "tolerance": self.tolerance,
            "min_distance": self.min_distance,
            "radii_source": self.radii_source,
            "tolerance_source": self.tolerance_source,
            "min_distance_source": self.min_distance_source,
            "metals_source": self.metals_source,
        }


def _section(data: dict[str, Any], key: str) -> dict[str, Any]:
    section = data.get(key)
    if not isinstance(section, dict):
        raise ReferenceDataError(f"{COVALENT_FILE} : section « {key} » absente")
    return section


def _text(section: dict[str, Any], key: str, where: str) -> str:
    value = section.get(key)
    if not value:
        raise ReferenceDataError(f"{COVALENT_FILE} : « {where}.{key} » absent")
    return str(value)


def _positive(section: dict[str, Any], key: str, where: str) -> float:
    value = section.get(key)
    if isinstance(value, bool) or not isinstance(value, int | float) or not value > 0:
        raise ReferenceDataError(
            f"{COVALENT_FILE} : « {where}.{key} » doit être un nombre positif, lu {value!r}"
        )
    return float(value)


def load_covalent_criteria() -> CovalentCriteria:
    """Charge et vérifie `covalent_geometry.yaml` ; chaque valeur doit avoir sa source."""
    data = read_data_file(COVALENT_FILE)
    radii_section = _section(data, "covalent_radii")
    detection = _section(data, "detection")
    metals_section = _section(data, "metals")
    values = radii_section.get("values")
    if not isinstance(values, dict) or not values:
        raise ReferenceDataError(f"{COVALENT_FILE} : « covalent_radii.values » absent ou vide")
    radii = {str(k).upper(): _positive(values, k, "covalent_radii.values") for k in values}
    metals_list = metals_section.get("elements")
    if not isinstance(metals_list, list) or not metals_list:
        raise ReferenceDataError(f"{COVALENT_FILE} : « metals.elements » absent ou vide")
    metals = frozenset(str(e).upper() for e in metals_list)
    both = sorted(metals & set(radii))
    if both:
        raise ReferenceDataError(
            f"{COVALENT_FILE} : élément(s) à la fois métal et rayon covalent : {both}"
        )
    tolerance = _positive(detection, "tolerance", "detection")
    min_distance = _positive(detection, "min_distance", "detection")
    criteria = CovalentCriteria(
        radii=radii,
        tolerance=tolerance,
        min_distance=min_distance,
        metals=metals,
        radii_source=_text(radii_section, "source", "covalent_radii"),
        tolerance_source=_text(detection, "tolerance_source", "detection"),
        min_distance_source=_text(detection, "min_distance_source", "detection"),
        metals_source=_text(metals_section, "source", "metals"),
    )
    smallest = 2 * min(radii.values()) + tolerance
    if not min_distance < smallest:
        raise ReferenceDataError(
            f"{COVALENT_FILE} : distance minimale {min_distance} Å supérieure à la plus courte "
            f"liaison possible ({smallest} Å) : aucune liaison ne serait reconnue"
        )
    return criteria
