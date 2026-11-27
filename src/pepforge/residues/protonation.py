"""Groupes ionisables, choix des états de protonation et charge nette théorique.

Les pKa, les états et leurs charges viennent de `data/ionizable_groups.yaml`,
où chaque valeur cite sa source. Rien n'est codé en dur ici, et rien n'est
supposé : un groupe sans pKa, un résidu ou une extrémité que le fichier ne
classe pas, un état qui n'appartient pas au groupe, tout cela lève une
exception au lieu d'être complété en silence.

Règle de choix d'un état, pour chaque groupe ionisable du peptide :

1. si un état est déjà fixé (configuration ou peptide d'entrée), il est
   **conservé** et marqué « imposé », après vérification qu'il appartient bien
   au groupe ; il peut être la forme minoritaire, ce que les contrôles signalent ;
2. sinon, la forme majoritaire au pH est retenue : forme acide (protonée) si
   pH < pKa, forme basique si pH > pKa ;
3. à pH = pKa exactement, aucune forme n'est majoritaire : l'état doit être
   imposé ;
4. si la forme majoritaire a plusieurs états (les deux tautomères neutres de
   l'histidine), le pKa ne permet pas de choisir : l'état doit être imposé.

Deux charges nettes sont calculées. La **charge des états retenus**, entière,
est celle que la topologie devra porter. La **charge théorique moyenne**
(`theoretical_net_charge`) suit Henderson–Hasselbalch : pour un groupe de pKa
donné, la fraction de forme acide vaut 1 / (1 + 10^(pH − pKa)).

Hors du périmètre de cette étape : le déplacement des pKa par l'environnement
(protéine repliée, complexe), qui demande un calcul sur une structure 3D.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from pepforge.core.peptide import Peptide
from pepforge.core.protonation import (
    C_TERMINUS_SITE,
    N_TERMINUS_SITE,
    GroupAssignment,
    IonizationForm,
    ProtonationResult,
)
from pepforge.exceptions import (
    ChemistryError,
    ConfigurationError,
    ReferenceDataError,
    UnsupportedFeatureError,
)
from pepforge.residues.library import read_data_file

DATA_FILE = "ionizable_groups.yaml"

PH_SCALE = (0.0, 14.0)
"""Plage de pH acceptée : l'échelle usuelle des solutions aqueuses diluées.

Ce n'est pas un paramètre scientifique mais une garde contre une faute de
saisie (un pH de 74 au lieu de 7,4). Les pKa de référence n'ont de sens que
dans l'eau, et bien en deçà de ces bornes.
"""


@dataclass(frozen=True)
class IonizationState:
    """Un état d'un groupe ionisable : clé interne, forme et charge formelle."""

    key: str
    form: IonizationForm
    charge: int


@dataclass(frozen=True)
class IonizableGroup:
    """Groupe ionisable porté par un résidu ou par une extrémité."""

    key: str
    label: str
    residue: str | None
    terminus: str | None
    pka: float | None
    """`None` : valeur manquante, tout peptide qui porte ce groupe est refusé."""

    pka_source: str
    states: tuple[IonizationState, ...]

    def states_of(self, form: IonizationForm) -> tuple[IonizationState, ...]:
        return tuple(s for s in self.states if s.form is form)

    def charge_of(self, form: IonizationForm) -> int:
        """Charge commune à tous les états d'une forme (vérifiée au chargement)."""
        return self.states_of(form)[0].charge

    def state(self, key: str) -> IonizationState | None:
        return next((s for s in self.states if s.key == key), None)

    def require_pka(self, where: str) -> float:
        if self.pka is None:
            raise ReferenceDataError(
                f"{where} : pKa de référence non renseigné pour le groupe {self.key} "
                f"({self.label}) dans {DATA_FILE} (source prévue : {self.pka_source}). "
                "La valeur doit être saisie avec sa source avant de traiter ce peptide."
            )
        return self.pka


@dataclass(frozen=True)
class IonizableGroupTable:
    """Contenu validé de `ionizable_groups.yaml`."""

    version: int
    data_set: str
    source: str
    sources: Mapping[str, str]
    """Identifiant → citation complète."""

    near_pka_window: float
    by_residue: Mapping[str, IonizableGroup]
    by_terminus: Mapping[str, IonizableGroup]
    non_ionizable_residues: frozenset[str]
    non_ionizable_termini: frozenset[str]

    @property
    def reference(self) -> str:
        """Référence courte du jeu de données, recopiée dans chaque résultat."""
        return f"{DATA_FILE} v{self.version} ({self.data_set})"

    def group_for_residue(self, code: str, where: str) -> IonizableGroup | None:
        group = self.by_residue.get(code)
        if group is None and code not in self.non_ionizable_residues:
            raise ReferenceDataError(
                f"{where} : le résidu {code} n'est classé ni ionisable ni non ionisable "
                f"dans {DATA_FILE}. Ses groupes ionisables doivent être décrits, avec "
                "leur source, avant de pouvoir choisir son état."
            )
        return group

    def group_for_terminus(self, key: str, where: str) -> IonizableGroup | None:
        group = self.by_terminus.get(key)
        if group is None and key not in self.non_ionizable_termini:
            raise ReferenceDataError(
                f"{where} : l'extrémité {key} n'est classée ni ionisable ni non ionisable "
                f"dans {DATA_FILE}."
            )
        return group


# --- Chargement et validation du fichier de données -------------------------


def _fail(message: str) -> ReferenceDataError:
    return ReferenceDataError(f"{DATA_FILE} invalide : {message}")


def _require(entry: Mapping[str, Any], key: str, kind: type, where: str) -> Any:
    if key not in entry:
        raise _fail(f"{where} : champ « {key} » absent")
    value = entry[key]
    if kind is float:
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise _fail(f"{where} : « {key} » doit être un nombre, reçu {value!r}")
        return float(value)
    if kind is int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise _fail(f"{where} : « {key} » doit être un entier, reçu {value!r}")
        return value
    if not isinstance(value, kind):
        raise _fail(f"{where} : « {key} » de type inattendu ({type(value).__name__})")
    return value


def _parse_states(entry: Mapping[str, Any], where: str) -> tuple[IonizationState, ...]:
    raw = _require(entry, "states", list, where)
    states: list[IonizationState] = []
    for item in raw:
        if not isinstance(item, dict):
            raise _fail(f"{where} : état mal formé {item!r}")
        key = _require(item, "key", str, where)
        try:
            form = IonizationForm(_require(item, "form", str, f"{where}, état {key}"))
        except ValueError:
            raise _fail(f"{where}, état {key} : forme inconnue {item['form']!r}") from None
        states.append(IonizationState(key, form, _require(item, "charge", int, where)))
    for form in IonizationForm:
        charges = {s.charge for s in states if s.form is form}
        if not charges:
            raise _fail(f"{where} : aucun état de forme {form.label}")
        if len(charges) > 1:
            raise _fail(f"{where} : les états de forme {form.label} n'ont pas tous la même charge")
    acid = next(s.charge for s in states if s.form is IonizationForm.ACID)
    base = next(s.charge for s in states if s.form is IonizationForm.BASE)
    if acid - base != 1:
        raise _fail(
            f"{where} : la forme protonée doit porter exactement une charge de plus que la "
            f"forme déprotonée (un proton), reçu {acid:+d} et {base:+d}"
        )
    return tuple(states)


def parse_ionizable_groups(data: Mapping[str, Any]) -> IonizableGroupTable:
    """Valide le contenu d'un fichier de groupes ionisables et le convertit.

    Séparée du chargement pour que les tests puissent valider des variantes du
    fichier. Toute incohérence lève `ReferenceDataError`.
    """
    version = _require(data, "version", int, "en-tête")
    data_set = _require(data, "data_set", str, "en-tête")
    source = _require(data, "source", str, "en-tête")
    raw_sources = _require(data, "sources", dict, "en-tête")
    sources: dict[str, str] = {}
    for key, value in raw_sources.items():
        if not isinstance(value, dict):
            raise _fail(f"source {key} mal formée")
        sources[str(key)] = _require(value, "citation", str, f"source {key}")
    window_entry = _require(data, "near_pka_window", dict, "en-tête")
    window = _require(window_entry, "value", float, "near_pka_window")
    _require(window_entry, "justification", str, "near_pka_window")
    if not window > 0:
        raise _fail("near_pka_window.value doit être strictement positive")

    by_residue: dict[str, IonizableGroup] = {}
    by_terminus: dict[str, IonizableGroup] = {}
    group_keys: set[str] = set()
    state_keys: set[str] = set()
    for entry in _require(data, "groups", list, "en-tête"):
        if not isinstance(entry, dict):
            raise _fail(f"groupe mal formé {entry!r}")
        key = _require(entry, "key", str, "groupe")
        where = f"groupe {key}"
        if key in group_keys:
            raise _fail(f"{where} déclaré deux fois")
        group_keys.add(key)
        residue = entry.get("residue")
        terminus = entry.get("terminus")
        if (residue is None) == (terminus is None):
            raise _fail(f"{where} : exactement un des champs « residue » et « terminus » attendu")
        if "pka" not in entry:
            raise _fail(f"{where} : champ « pka » absent (écrire null si la valeur manque)")
        pka = None if entry["pka"] is None else _require(entry, "pka", float, where)
        pka_source = _require(entry, "pka_source", str, where)
        if pka_source not in sources:
            raise _fail(f"{where} : source {pka_source!r} non déclarée dans « sources »")
        states = _parse_states(entry, where)
        for state in states:
            if state.key in state_keys:
                raise _fail(f"{where} : clé d'état {state.key} déjà utilisée")
            state_keys.add(state.key)
        group = IonizableGroup(
            key=key,
            label=_require(entry, "label", str, where),
            residue=residue,
            terminus=terminus,
            pka=pka,
            pka_source=pka_source,
            states=states,
        )
        target = by_residue if residue is not None else by_terminus
        site = str(residue if residue is not None else terminus)
        if site in target:
            raise _fail(f"{where} : {site} porte déjà le groupe {target[site].key}")
        target[site] = group

    non_res = frozenset(str(c) for c in _require(data, "non_ionizable_residues", list, "en-tête"))
    non_ter = frozenset(str(t) for t in _require(data, "non_ionizable_termini", list, "en-tête"))
    overlaps = (
        ("résidu", non_res & set(by_residue)),
        ("extrémité", non_ter & set(by_terminus)),
    )
    for name, both in overlaps:
        if both:
            raise _fail(f"{name}(s) classé(s) ionisable(s) et non ionisable(s) : {sorted(both)}")
    return IonizableGroupTable(
        version=version,
        data_set=data_set,
        source=source,
        sources=sources,
        near_pka_window=window,
        by_residue=by_residue,
        by_terminus=by_terminus,
        non_ionizable_residues=non_res,
        non_ionizable_termini=non_ter,
    )


def load_ionizable_groups() -> IonizableGroupTable:
    """Charge et valide le fichier de groupes ionisables livré avec PepForge."""
    return parse_ionizable_groups(read_data_file(DATA_FILE))


# --- Calculs ------------------------------------------------------------------


def check_ph(ph: object) -> float:
    """Vérifie qu'un pH est un nombre fini dans l'échelle `PH_SCALE`."""
    if isinstance(ph, bool) or not isinstance(ph, int | float) or not math.isfinite(ph):
        raise ConfigurationError(f"pH invalide : {ph!r} (nombre attendu)")
    low, high = PH_SCALE
    if not low <= ph <= high:
        raise ConfigurationError(f"pH {ph} hors de l'échelle acceptée ({low} à {high})")
    return float(ph)


def acid_fraction(pka: float, ph: float) -> float:
    """Fraction de la forme acide (protonée) à l'équilibre, Henderson–Hasselbalch."""
    return 1.0 / (1.0 + math.pow(10.0, ph - pka))


def major_form(pka: float, ph: float) -> IonizationForm | None:
    """Forme majoritaire au pH, `None` à pH = pKa exactement."""
    if ph < pka:
        return IonizationForm.ACID
    if ph > pka:
        return IonizationForm.BASE
    return None


@dataclass(frozen=True)
class _Site:
    site: str
    residue: str
    group: IonizableGroup
    current_state: str | None

    @property
    def where(self) -> str:
        if self.site == N_TERMINUS_SITE:
            return "extrémité N"
        if self.site == C_TERMINUS_SITE:
            return "extrémité C"
        return f"résidu {self.residue} {self.site}"


def _refuse_unsupported(peptide: Peptide) -> None:
    unsupported: list[str] = []
    if peptide.is_cyclic:
        unsupported.append("cyclisation tête-queue (les groupes terminaux sont consommés)")
    unsupported += [
        f"{m.describe()} (elle peut ajouter ou retirer un groupe ionisable)"
        for m in peptide.modifications
    ]
    if unsupported:
        raise UnsupportedFeatureError(
            "La protonation ne sait pas encore traiter : " + "; ".join(unsupported)
        )


def _sites(peptide: Peptide, table: IonizableGroupTable) -> list[_Site]:
    """Groupes ionisables du peptide, dans l'ordre N, résidus, C.

    Une cystéine engagée dans un pont disulfure n'a plus de thiol : elle n'est
    pas ionisable. Un état fixé sur un site sans groupe ionisable est une
    incohérence chimique, signalée et jamais effacée.
    """
    _refuse_unsupported(peptide)
    bridged = {p for pair in peptide.disulfides for p in pair}
    sites: list[_Site] = []

    def add_terminus(side: str, key: str, state: str | None, where: str) -> None:
        group = table.group_for_terminus(key, where)
        if group is not None:
            sites.append(_Site(side, key, group, state))
        elif state is not None:
            raise ChemistryError(
                f"{where} ({key}) : état de protonation {state!r} fixé sur une extrémité "
                "sans groupe ionisable"
            )

    n_terminus, c_terminus = peptide.n_terminus, peptide.c_terminus
    add_terminus(N_TERMINUS_SITE, n_terminus.key, n_terminus.protonation, "extrémité N")
    for position, residue in enumerate(peptide.residues, start=1):
        where = f"résidu {residue.code} {position}"
        if position in bridged:
            if residue.protonation is not None:
                raise ChemistryError(
                    f"{where} : état {residue.protonation!r} fixé sur un résidu engagé dans "
                    "un pont disulfure, qui n'a plus de groupe ionisable"
                )
            continue
        group = table.group_for_residue(residue.code, where)
        if group is not None:
            sites.append(_Site(str(position), residue.code, group, residue.protonation))
        elif residue.protonation is not None:
            raise ChemistryError(
                f"{where} : état {residue.protonation!r} fixé sur un résidu sans groupe ionisable"
            )
    add_terminus(C_TERMINUS_SITE, c_terminus.key, c_terminus.protonation, "extrémité C")
    return sites


def theoretical_net_charge(
    peptide: Peptide, ph: float, table: IonizableGroupTable | None = None
) -> float:
    """Charge nette moyenne au pH, somme des charges moyennes de chaque groupe.

    Les états éventuellement fixés sur le peptide n'entrent pas dans ce calcul :
    c'est la charge d'équilibre au pH, indépendante de tout choix.
    """
    ph = check_ph(ph)
    table = table or load_ionizable_groups()
    total = 0.0
    for site in _sites(peptide, table):
        group = site.group
        fraction = acid_fraction(group.require_pka(site.where), ph)
        acid = group.charge_of(IonizationForm.ACID)
        base = group.charge_of(IonizationForm.BASE)
        total += base + (acid - base) * fraction
    return total


def _choose(site: _Site, ph: float, major: IonizationForm | None) -> tuple[IonizationState, bool]:
    group = site.group
    if site.current_state is not None:
        state = group.state(site.current_state)
        if state is None:
            known = ", ".join(s.key for s in group.states)
            raise ChemistryError(
                f"{site.where} : état {site.current_state!r} inconnu pour le groupe "
                f"{group.key} (états possibles : {known})"
            )
        return state, True
    if major is None:
        raise ChemistryError(
            f"{site.where} : pH {ph} égal au pKa ({group.pka}) du groupe {group.key}, aucune "
            "forme n'est majoritaire. L'état doit être imposé dans la configuration."
        )
    candidates = group.states_of(major)
    if len(candidates) > 1:
        keys = " ou ".join(s.key for s in candidates)
        raise ChemistryError(
            f"{site.where} : forme {major.label} majoritaire à pH {ph}, mais le pKa ne "
            f"départage pas ses états ({keys}). L'état doit être imposé dans la configuration."
        )
    return candidates[0], False


def assign_protonation(
    peptide: Peptide, ph: float, table: IonizableGroupTable | None = None
) -> ProtonationResult:
    """Choisit l'état de chaque groupe ionisable et calcule les charges nettes.

    Renvoie un `ProtonationResult` dont le peptide porte les états retenus
    (`Residue.protonation`, `Terminus.protonation`), avec la trace de chaque
    choix. Le peptide d'entrée n'est pas modifié ; la journalisation revient à
    l'appelant (`pepforge.pipeline.ProtonationStage`).
    """
    ph = check_ph(ph)
    table = table or load_ionizable_groups()
    assignments: list[GroupAssignment] = []
    result_peptide = peptide
    for site in _sites(peptide, table):
        pka = site.group.require_pka(site.where)
        major = major_form(pka, ph)
        state, imposed = _choose(site, ph, major)
        assignments.append(
            GroupAssignment(
                site=site.site,
                residue=site.residue,
                group=site.group.key,
                state=state.key,
                form=state.form,
                major_form=major,
                charge=state.charge,
                pka=pka,
                pka_source=site.group.pka_source,
                imposed=imposed,
            )
        )
        if site.site == N_TERMINUS_SITE:
            result_peptide = result_peptide.with_n_terminus_protonation(state.key)
        elif site.site == C_TERMINUS_SITE:
            result_peptide = result_peptide.with_c_terminus_protonation(state.key)
        else:
            result_peptide = result_peptide.with_protonation(int(site.site), state.key)
    return ProtonationResult(
        peptide=result_peptide,
        ph=ph,
        assignments=tuple(assignments),
        net_charge=sum(a.charge for a in assignments),
        theoretical_net_charge=theoretical_net_charge(peptide, ph, table),
        near_pka_window=table.near_pka_window,
        data_source=table.reference,
    )
