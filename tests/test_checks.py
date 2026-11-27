"""Tests des contrôles qualité et de l'exécution contrôlée des étapes."""

from __future__ import annotations

from pathlib import Path

import pytest

from pepforge.checks import PEPTIDE_CHECKS, Check, CheckResult, Status, run_checks
from pepforge.config import PrepConfig
from pepforge.core import Link, LinkKind, Modification, Peptide, Residue
from pepforge.exceptions import QualityCheckError
from pepforge.io import SequenceNotation, parse_sequence
from pepforge.pipeline import SequenceInputStage, StageContext, run_stage


def _context(tmp_path: Path) -> StageContext:
    return StageContext(config=PrepConfig(), workdir=tmp_path)


def test_declared_termini_pass() -> None:
    results = run_checks(PEPTIDE_CHECKS, parse_sequence("Ac-RGDS-NH2", notation="one_letter"))
    assert all(r.status is Status.PASS for r in results)


def test_default_termini_warn_but_do_not_block(tmp_path: Path) -> None:
    context = _context(tmp_path)
    peptide = run_stage(
        SequenceInputStage(SequenceNotation.ONE_LETTER), "RGDS", context, PEPTIDE_CHECKS
    )
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
        run_stage(
            SequenceInputStage(SequenceNotation.ONE_LETTER), "Ac-RGDS-NH2", context, [always_fail]
        )
    # Les résultats restent disponibles pour le rapport, même en cas d'échec.
    assert context.check_results["input.sequence"][0].status is Status.FAIL


def test_stage_journals_the_parse(tmp_path: Path) -> None:
    context = _context(tmp_path)
    run_stage(SequenceInputStage(SequenceNotation.ONE_LETTER), "Ac-RGDS-NH2", context)
    actions = [e.action for e in context.journal.entries]
    assert actions == ["parse_sequence", "stage_completed"]


# --- Liaisons et modifications déclarées (étape 1) ------------------------------


def _check(check_id: str, peptide: Peptide) -> CheckResult:
    results = {r.check_id: r for r in run_checks(PEPTIDE_CHECKS, peptide)}
    return results[check_id]


def test_disulfide_between_two_cysteines_passes() -> None:
    peptide = parse_sequence("Ac-CRGDC-NH2", notation="one_letter").with_link(
        Link(LinkKind.DISULFIDE, 1, 5)
    )
    result = _check("peptide.links_chemistry", peptide)
    assert result.status is Status.PASS
    assert "pont disulfure 1–5" in result.message


def test_disulfide_on_a_residue_without_thiol_fails() -> None:
    peptide = parse_sequence("Ac-CRGDS-NH2", notation="one_letter").with_link(
        Link(LinkKind.DISULFIDE, 1, 5)
    )
    result = _check("peptide.links_chemistry", peptide)
    assert result.status is Status.FAIL
    assert "SER en position 5 ne porte pas de thiol" in result.message


def test_disulfide_on_a_non_standard_residue_is_unverifiable_hence_fails() -> None:
    xcy = Residue(code="XCY", standard=False)
    peptide = parse_sequence("Ac-CRGDA-NH2", notation="one_letter").mutate(5, xcy)
    peptide = peptide.with_link(Link(LinkKind.DISULFIDE, 1, 5))
    result = _check("peptide.links_chemistry", peptide)
    assert result.status is Status.FAIL
    assert "non vérifiable" in result.message


def test_head_to_tail_with_free_ends_passes() -> None:
    cyclic = parse_sequence("H-RGDS-OH", notation="one_letter").with_link(
        Link(LinkKind.HEAD_TO_TAIL, 1, 4)
    )
    assert _check("peptide.links_chemistry", cyclic).status is Status.PASS


def test_head_to_tail_with_capped_ends_fails() -> None:
    cyclic = parse_sequence("Ac-RGDS-NH2", notation="one_letter").with_link(
        Link(LinkKind.HEAD_TO_TAIL, 1, 4)
    )
    result = _check("peptide.links_chemistry", cyclic)
    assert result.status is Status.FAIL
    assert "extrémité N" in result.message
    assert "extrémité C" in result.message


def test_head_to_tail_on_a_non_standard_end_is_unverifiable_hence_fails() -> None:
    cyclic = parse_sequence("H-RGDS-OH", notation="one_letter").mutate(
        4, Residue(code="XAA", standard=False)
    )
    cyclic = cyclic.with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 4))
    result = _check("peptide.links_chemistry", cyclic)
    assert result.status is Status.FAIL
    assert "non vérifiable" in result.message


def test_declared_modification_is_unverifiable_hence_fails() -> None:
    peptide = parse_sequence("Ac-RGDS-NH2", notation="one_letter").with_modification(
        Modification("phosphorylation", 4)
    )
    result = _check("peptide.modifications_verifiable", peptide)
    assert result.status is Status.FAIL
    assert "phosphorylation en position 4" in result.message


def test_a_problem_on_a_link_stops_the_input_stage(tmp_path: Path) -> None:
    class _CappedCyclicStage(SequenceInputStage):
        def run(self, data: str, context: StageContext) -> Peptide:
            peptide = super().run(data, context)
            return peptide.with_link(Link(LinkKind.HEAD_TO_TAIL, 1, len(peptide)))

    with pytest.raises(QualityCheckError, match=r"peptide\.links_chemistry"):
        run_stage(
            _CappedCyclicStage(SequenceNotation.ONE_LETTER),
            "Ac-RGDS-NH2",
            _context(tmp_path),
            PEPTIDE_CHECKS,
        )


@pytest.mark.parametrize(
    ("notation", "length"),
    [(SequenceNotation.ONE_LETTER, 3), (SequenceNotation.THREE_LETTER, 1)],
)
def test_input_stage_reads_in_the_declared_notation_and_journals_it(
    tmp_path: Path, notation: SequenceNotation, length: int
) -> None:
    """`Ac-ALA-NHMe` : Ala-Leu-Ala en code à une lettre, l'alanine seule en code à trois lettres."""
    context = _context(tmp_path)
    peptide = run_stage(SequenceInputStage(notation), "Ac-ALA-NHMe", context, PEPTIDE_CHECKS)
    assert len(peptide) == length
    entry = context.journal.entries[0]
    assert entry.action == "parse_sequence"
    assert entry.parameters["notation"] == notation.value
    assert type(entry.parameters["notation"]) is str  # sérialisable par yaml.safe_dump
    statuses = {r.check_id: r.status for r in context.check_results["input.sequence"]}
    assert all(status is Status.PASS for status in statuses.values())


def test_link_checks_point_to_the_checklist_item() -> None:
    results = {
        r.check_id: r
        for r in run_checks(PEPTIDE_CHECKS, parse_sequence("Ac-RGDS-NH2", notation="one_letter"))
    }
    for check_id in ("peptide.links_chemistry", "peptide.modifications_verifiable"):
        assert results[check_id].checklist_item == "structure.links_and_modifications"
