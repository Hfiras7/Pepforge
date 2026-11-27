"""Contrôles appliqués au `ProtonationResult` à la sortie de l'étape de protonation.

Ils ne refont pas le choix des états : ils vérifient que le résultat est
cohérent avec lui-même (chaque état décidé est bien porté par le peptide, la
charge nette est la somme des charges des états, chaque pKa a sa source), et
signalent ce qui doit être relu (pH proche d'un pKa, état imposé minoritaire).
"""

from __future__ import annotations

from pepforge.checks.base import Check, CheckResult, Status
from pepforge.core.protonation import (
    C_TERMINUS_SITE,
    N_TERMINUS_SITE,
    GroupAssignment,
    ProtonationResult,
)

_HYDROGENS_ITEM = "structure.hydrogens_ph"
_CHARGE_ITEM = "parametrisation.net_charge"


def _carried_states(result: ProtonationResult) -> dict[str, str | None]:
    peptide = result.peptide
    states: dict[str, str | None] = {
        N_TERMINUS_SITE: peptide.n_terminus.protonation,
        C_TERMINUS_SITE: peptide.c_terminus.protonation,
    }
    for position, residue in enumerate(peptide.residues, start=1):
        states[str(position)] = residue.protonation
    return states


def _states_recorded(result: ProtonationResult) -> CheckResult:
    carried = _carried_states(result)
    decided = {a.site: a.state for a in result.assignments}
    problems = [
        f"site {site} : {decided.get(site)!r} décidé, {carried.get(site)!r} porté par le peptide"
        for site in sorted(set(carried) | set(decided))
        if carried.get(site) != decided.get(site)
    ]
    if problems:
        return CheckResult(
            "protonation.states_recorded",
            Status.FAIL,
            "États décidés et états portés par le peptide divergent : " + "; ".join(problems),
            _HYDROGENS_ITEM,
        )
    return CheckResult(
        "protonation.states_recorded",
        Status.PASS,
        f"{len(decided)} groupe(s) ionisable(s), chacun avec son état",
        _HYDROGENS_ITEM,
    )


def _net_charge(result: ProtonationResult) -> CheckResult:
    total = sum(a.charge for a in result.assignments)
    if isinstance(result.net_charge, bool) or not isinstance(result.net_charge, int):
        return CheckResult(
            "protonation.net_charge",
            Status.FAIL,
            f"Charge nette non entière : {result.net_charge!r}",
            _CHARGE_ITEM,
        )
    if total != result.net_charge:
        return CheckResult(
            "protonation.net_charge",
            Status.FAIL,
            f"Charge nette annoncée {result.net_charge:+d}, somme des états {total:+d}",
            _CHARGE_ITEM,
        )
    return CheckResult(
        "protonation.net_charge",
        Status.PASS,
        f"Charge nette {result.net_charge:+d} à pH {result.ph} "
        f"(charge moyenne au pH : {result.theoretical_net_charge:+.3f})",
        _CHARGE_ITEM,
        {"net_charge": str(result.net_charge), "ph": str(result.ph)},
    )


def _sources(result: ProtonationResult) -> CheckResult:
    unsourced = [a.describe() for a in result.assignments if not a.pka_source.strip()]
    if not result.data_source.strip():
        unsourced.append("jeu de données sans référence")
    if unsourced:
        return CheckResult(
            "protonation.sources",
            Status.FAIL,
            "pKa sans source : " + "; ".join(unsourced),
            _HYDROGENS_ITEM,
        )
    return CheckResult(
        "protonation.sources",
        Status.PASS,
        f"Chaque pKa cite sa source ({result.data_source})",
        _HYDROGENS_ITEM,
    )


def _deduced_states_follow_ph(result: ProtonationResult) -> CheckResult:
    wrong = [
        f"{a.describe()}, forme {a.form.label} alors que la forme "
        f"{a.major_form.label if a.major_form else 'majoritaire (aucune)'} domine"
        for a in result.assignments
        if not a.imposed and a.form is not a.major_form
    ]
    if wrong:
        return CheckResult(
            "protonation.deduced_states_follow_ph",
            Status.FAIL,
            "État non imposé différent de la forme majoritaire : " + "; ".join(wrong),
            _HYDROGENS_ITEM,
        )
    return CheckResult(
        "protonation.deduced_states_follow_ph",
        Status.PASS,
        "Chaque état non imposé est la forme majoritaire au pH",
        _HYDROGENS_ITEM,
    )


def _margin(assignment: GroupAssignment, ph: float) -> float:
    return abs(ph - assignment.pka)


def _near_pka(result: ProtonationResult) -> CheckResult:
    close = [
        f"{a.describe()}, pKa {a.pka} (écart {_margin(a, result.ph):.2f})"
        for a in result.assignments
        if _margin(a, result.ph) < result.near_pka_window
    ]
    if close:
        return CheckResult(
            "protonation.near_pka",
            Status.WARN,
            f"pH {result.ph} à moins de {result.near_pka_window} unité du pKa, état à relire "
            "(les deux formes coexistent en proportion notable) : " + "; ".join(close),
            _HYDROGENS_ITEM,
        )
    return CheckResult(
        "protonation.near_pka",
        Status.PASS,
        f"Aucun pKa à moins de {result.near_pka_window} unité du pH {result.ph}",
        _HYDROGENS_ITEM,
    )


def _imposed_minority(result: ProtonationResult) -> CheckResult:
    minority = [
        a.describe() for a in result.assignments if a.imposed and a.form is not a.major_form
    ]
    if minority:
        return CheckResult(
            "protonation.imposed_minority",
            Status.WARN,
            f"État imposé qui n'est pas la forme majoritaire à pH {result.ph} : "
            + "; ".join(minority)
            + ". Vérifier que ce choix est voulu et justifié dans la configuration.",
            _HYDROGENS_ITEM,
        )
    imposed = sum(a.imposed for a in result.assignments)
    return CheckResult(
        "protonation.imposed_minority",
        Status.PASS,
        f"{imposed} état(s) imposé(s), tous de forme majoritaire",
        _HYDROGENS_ITEM,
    )


PROTONATION_CHECKS: tuple[Check[ProtonationResult], ...] = (
    Check(
        "protonation.states_recorded",
        "chaque état décidé est porté par le peptide, et aucun autre",
        _states_recorded,
        _HYDROGENS_ITEM,
    ),
    Check(
        "protonation.net_charge",
        "charge nette entière, égale à la somme des charges des états",
        _net_charge,
        _CHARGE_ITEM,
    ),
    Check("protonation.sources", "chaque pKa cite sa source", _sources, _HYDROGENS_ITEM),
    Check(
        "protonation.deduced_states_follow_ph",
        "un état non imposé est la forme majoritaire au pH",
        _deduced_states_follow_ph,
        _HYDROGENS_ITEM,
    ),
    Check(
        "protonation.near_pka",
        "pH proche d'un pKa : état à relire",
        _near_pka,
        _HYDROGENS_ITEM,
    ),
    Check(
        "protonation.imposed_minority",
        "état imposé de forme minoritaire",
        _imposed_minority,
        _HYDROGENS_ITEM,
    ),
)
