"""Contrôles appliqués au `Peptide` à la sortie de l'étape d'entrée."""

from __future__ import annotations

from pepforge.checks.base import Check, CheckResult, Status
from pepforge.core.peptide import Peptide


def _has_residues(peptide: Peptide) -> CheckResult:
    if len(peptide) == 0:
        return CheckResult("peptide.non_empty", Status.FAIL, "Le peptide ne contient aucun résidu")
    return CheckResult(
        "peptide.non_empty", Status.PASS, f"{len(peptide)} résidu(s)", "structure.residue_count"
    )


def _termini_declared(peptide: Peptide) -> CheckResult:
    defaulted = [
        side
        for side, key in (("N", "n_terminus_source"), ("C", "c_terminus_source"))
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


PEPTIDE_CHECKS: tuple[Check[Peptide], ...] = (
    Check("peptide.non_empty", "Le peptide contient au moins un résidu", _has_residues),
    Check(
        "peptide.termini_declared",
        "Les deux extrémités sont déclarées explicitement",
        _termini_declared,
        "structure.termini",
    ),
)
"""Contrôles exécutés après la lecture d'une séquence."""
