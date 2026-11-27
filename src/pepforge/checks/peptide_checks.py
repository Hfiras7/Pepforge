"""Contrôles appliqués au `Peptide` à la sortie de l'étape d'entrée."""

from __future__ import annotations

from typing import assert_never

from pepforge.checks.base import Check, CheckResult, Status
from pepforge.core.peptide import (
    C_TERMINUS_DEFAULTED,
    N_TERMINUS_DEFAULTED,
    Link,
    LinkKind,
    Peptide,
)

_THIOL_RESIDUES = frozenset({"CYS"})
"""Codes des résidus dont la chaîne latérale porte un thiol : la cystéine.

Liste volontairement fermée. Un pont disulfure déclaré sur un autre résidu
échoue : il est impossible sur un résidu standard sans thiol, et non vérifiable
sur un résidu non standard, dont la bibliothèque ne décrit pas encore les
groupes fonctionnels.
"""

# Une cyclisation tête-queue consomme l'amine α libre du premier résidu et le
# carboxyle α libre du dernier : seules ces deux extrémités lui sont compatibles.
_HEAD_TO_TAIL_N_TERMINUS = "free_amine"
_HEAD_TO_TAIL_C_TERMINUS = "free_acid"


def _has_residues(peptide: Peptide) -> CheckResult:
    if len(peptide) == 0:
        return CheckResult("peptide.non_empty", Status.FAIL, "Le peptide ne contient aucun résidu")
    return CheckResult(
        "peptide.non_empty", Status.PASS, f"{len(peptide)} résidu(s)", "structure.residue_count"
    )


def _termini_declared(peptide: Peptide) -> CheckResult:
    defaulted = [
        side
        for side, key in (("N", N_TERMINUS_DEFAULTED), ("C", C_TERMINUS_DEFAULTED))
        if key in peptide.metadata
    ]
    if defaulted:
        return CheckResult(
            "peptide.termini_declared",
            Status.WARN,
            "Extrémité(s) non déclarée(s), valeur par défaut appliquée : " + ", ".join(defaulted),
            "structure.termini",
        )
    return CheckResult(
        "peptide.termini_declared",
        Status.PASS,
        f"N : {peptide.n_terminus.label} ; C : {peptide.c_terminus.label}",
        "structure.termini",
    )


_LINKS_ITEM = "structure.links_and_modifications"
"""Item de la checklist 04 sur les liaisons hors squelette et les modifications."""


def _disulfide_problems(peptide: Peptide, link: Link) -> list[str]:
    problems: list[str] = []
    for position in link.positions:
        residue = peptide.residue_at(position)
        if residue.code in _THIOL_RESIDUES:
            continue
        if residue.standard:
            problems.append(
                f"{link.describe()} : {residue.code} en position {position} ne porte pas de thiol"
            )
        else:
            problems.append(
                f"{link.describe()} : non vérifiable, la bibliothèque ne dit pas si "
                f"{residue.code} (position {position}) porte un thiol"
            )
    return problems


def _head_to_tail_problems(peptide: Peptide, link: Link) -> list[str]:
    problems: list[str] = []
    if peptide.n_terminus.key != _HEAD_TO_TAIL_N_TERMINUS:
        problems.append(
            f"{link.describe()} : l'amine α du premier résidu est consommée par la "
            f"cyclisation, l'extrémité N ne peut pas être « {peptide.n_terminus.label} »"
        )
    if peptide.c_terminus.key != _HEAD_TO_TAIL_C_TERMINUS:
        problems.append(
            f"{link.describe()} : le carboxyle α du dernier résidu est consommé par la "
            f"cyclisation, l'extrémité C ne peut pas être « {peptide.c_terminus.label} »"
        )
    for position in link.positions:
        residue = peptide.residue_at(position)
        if not residue.standard:
            problems.append(
                f"{link.describe()} : non vérifiable, la bibliothèque ne décrit pas les "
                f"groupes terminaux de {residue.code} (position {position})"
            )
    return problems


def _links_chemistry(peptide: Peptide) -> CheckResult:
    problems: list[str] = []
    for link in peptide.links:
        kind = link.kind
        if kind is LinkKind.DISULFIDE:
            problems += _disulfide_problems(peptide, link)
        elif kind is LinkKind.HEAD_TO_TAIL:
            problems += _head_to_tail_problems(peptide, link)
        else:
            assert_never(kind)
    if problems:
        return CheckResult("peptide.links_chemistry", Status.FAIL, "; ".join(problems), _LINKS_ITEM)
    if not peptide.links:
        return CheckResult(
            "peptide.links_chemistry",
            Status.PASS,
            "Aucune liaison hors du squelette linéaire",
            _LINKS_ITEM,
        )
    return CheckResult(
        "peptide.links_chemistry",
        Status.PASS,
        "Liaisons cohérentes : " + ", ".join(link.describe() for link in peptide.links),
        _LINKS_ITEM,
    )


def _modifications_verifiable(peptide: Peptide) -> CheckResult:
    if not peptide.modifications:
        return CheckResult(
            "peptide.modifications_verifiable",
            Status.PASS,
            "Aucune modification de chaîne latérale déclarée",
            _LINKS_ITEM,
        )
    declared = ", ".join(m.describe() for m in peptide.modifications)
    return CheckResult(
        "peptide.modifications_verifiable",
        Status.FAIL,
        f"Non vérifiable : aucune bibliothèque de modifications n'est encore chargée ({declared})",
        _LINKS_ITEM,
    )


PEPTIDE_CHECKS: tuple[Check[Peptide], ...] = (
    Check("peptide.non_empty", "Le peptide contient au moins un résidu", _has_residues),
    Check(
        "peptide.termini_declared",
        "Les deux extrémités sont déclarées explicitement",
        _termini_declared,
        "structure.termini",
    ),
    Check(
        "peptide.links_chemistry",
        "Les liaisons déclarées sont chimiquement possibles et vérifiables",
        _links_chemistry,
        _LINKS_ITEM,
    ),
    Check(
        "peptide.modifications_verifiable",
        "Les modifications de chaîne latérale déclarées sont vérifiables",
        _modifications_verifiable,
        _LINKS_ITEM,
    ),
)
"""Contrôles exécutés après la lecture d'une séquence."""
