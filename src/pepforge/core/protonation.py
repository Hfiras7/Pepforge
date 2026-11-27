"""Résultat de l'étape de protonation : état choisi pour chaque groupe ionisable.

Ces types ne contiennent que des données décrites ailleurs : les pKa, les
charges des états et leurs sources viennent de
`pepforge/residues/data/ionizable_groups.yaml`, chargé par
`pepforge.residues.protonation`. Ils vivent dans `core` pour que `checks`
puisse les contrôler sans connaître la bibliothèque de résidus (matrice de
dépendances, `docs/ARCHITECTURE.md`, section 3). Aucun nom de champ de force
n'y figure : les états sont des clés internes neutres (`his_epsilon`,
`nterm_protonated`...).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from pepforge.core.peptide import Peptide

N_TERMINUS_SITE = "N"
"""Valeur de `GroupAssignment.site` pour le groupe de l'extrémité N."""

C_TERMINUS_SITE = "C"
"""Valeur de `GroupAssignment.site` pour le groupe de l'extrémité C."""


class IonizationForm(StrEnum):
    """Les deux formes d'un couple acide–base."""

    ACID = "acid"
    """Forme protonée (porte le proton échangeable)."""

    BASE = "base"
    """Forme déprotonée."""

    @property
    def label(self) -> str:
        return "protonée" if self is IonizationForm.ACID else "déprotonée"


@dataclass(frozen=True)
class GroupAssignment:
    """État retenu pour un groupe ionisable du peptide, avec sa justification."""

    site: str
    """`"N"`, `"C"`, ou la position du résidu (à partir de 1) écrite en chiffres."""

    residue: str
    """Code du résidu porteur, ou clé de l'extrémité (`free_amine`...)."""

    group: str
    """Clé du groupe ionisable dans le fichier de données (`asp_side_chain`...)."""

    state: str
    """Clé interne de l'état retenu (`asp_deprotonated`...)."""

    form: IonizationForm
    """Forme de l'état retenu."""

    major_form: IonizationForm | None
    """Forme majoritaire au pH de l'étude ; `None` si pH = pKa exactement."""

    charge: int
    """Charge formelle de l'état retenu, en charge élémentaire."""

    pka: float
    pka_source: str
    """Identifiant de la source du pKa, défini dans le fichier de données."""

    imposed: bool
    """Vrai si l'état a été imposé (configuration ou peptide d'entrée), pas déduit du pH."""

    @property
    def is_residue(self) -> bool:
        return self.site not in (N_TERMINUS_SITE, C_TERMINUS_SITE)

    @property
    def position(self) -> int | None:
        """Position du résidu porteur, `None` pour un groupe terminal."""
        return int(self.site) if self.is_residue else None

    def describe(self) -> str:
        where = {N_TERMINUS_SITE: "extrémité N", C_TERMINUS_SITE: "extrémité C"}.get(
            self.site, f"{self.residue} {self.site}"
        )
        return f"{where} : {self.state} ({self.charge:+d})"


@dataclass(frozen=True)
class ProtonationResult:
    """Peptide protoné, charges nettes, et trace de chaque choix.

    `net_charge` est la somme des charges formelles des états retenus : c'est la
    charge entière que la topologie devra porter, référence des contrôles de
    charge ultérieurs. `theoretical_net_charge` est la charge moyenne au pH
    (Henderson–Hasselbalch), qui n'est pas entière : l'écart entre les deux dit
    à quel point le choix d'un état unique par groupe est une approximation.
    """

    peptide: Peptide
    ph: float
    assignments: tuple[GroupAssignment, ...]
    net_charge: int
    theoretical_net_charge: float
    near_pka_window: float
    """Fenêtre |pH − pKa| sous laquelle un état est signalé comme à relire."""

    data_source: str
    """Référence du jeu de pKa utilisé (fichier, version, identifiant)."""
