"""Tests des contrôles qualité et de l'exécution contrôlée des étapes."""

from __future__ import annotations

from pathlib import Path

import pytest

from pepforge.checks import PEPTIDE_CHECKS, Check, CheckResult, Status, run_checks
from pepforge.config import PrepConfig
from pepforge.core import Peptide
from pepforge.exceptions import QualityCheckError
from pepforge.io import parse_sequence
from pepforge.pipeline import SequenceInputStage, StageContext, run_stage


def _context(tmp_path: Path) -> StageContext:
    return StageContext(config=PrepConfig(), workdir=tmp_path)


def test_declared_termini_pass() -> None:
    results = run_checks(PEPTIDE_CHECKS, parse_sequence("Ac-RGDS-NH2"))
    assert all(r.status is Status.PASS for r in results)


def test_default_termini_warn_but_do_not_block(tmp_path: Path) -> None:
    context = _context(tmp_path)
    peptide = run_stage(SequenceInputStage(), "RGDS", context, PEPTIDE_CHECKS)
    assert len(peptide) == 4
    statuses = {r.check_id: r.status for r in context.check_results["input.sequence"]}
    assert statuses["peptide.termini_declared"] is Status.WARN
    assert context.journal.entries[-1].action == "stage_completed"


def test_failing_check_stops_the_pipeline(tmp_path: Path) -> None:
    always_fail: Check[Peptide] = Check(
        "test.fail", "échec forcé", lambda _: CheckResult("test.fail", Status.FAIL, "forcé")
    )
    context = _context(tmp_path)
    with pytest.raises(QualityCheckError, match=r"test\.fail"):
        run_stage(SequenceInputStage(), "Ac-RGDS-NH2", context, [always_fail])
    # Les résultats restent disponibles pour le rapport, même en cas d'échec.
    assert context.check_results["input.sequence"][0].status is Status.FAIL


def test_stage_journals_the_parse(tmp_path: Path) -> None:
    context = _context(tmp_path)
    run_stage(SequenceInputStage(), "Ac-RGDS-NH2", context)
    actions = [e.action for e in context.journal.entries]
    assert actions == ["parse_sequence", "stage_completed"]
