"""Tests du choix du niveau de validation et de l'audit de ce qui a été exécuté.

L'exigence vérifiée ici : demander un niveau et ne pas l'obtenir est une
erreur rapportée, jamais un silence.
"""

from __future__ import annotations

from pepforge.checks import CheckResult, Status, audit, requested, unavailable
from pepforge.checks.base import failures
from pepforge.checks.levels import CONTROLS
from pepforge.config import ValidationLevel


def _produced(*check_ids: str) -> list[CheckResult]:
    return [
        CheckResult(check_id=cid, status=Status.PASS, message="contrôle satisfait")
        for cid in check_ids
    ]


def test_control_ids_are_unique() -> None:
    ids = [c.control_id for c in CONTROLS]
    assert len(ids) == len(set(ids))


def test_levels_are_cumulative() -> None:
    level0 = requested(ValidationLevel.FILE_CONFORMITY)
    level3 = requested(ValidationLevel.ENERGETIC_COHERENCE)
    assert set(level0) <= set(level3)
    assert len(level3) > len(level0), "un niveau supérieur ajoute des contrôles"


def test_level_zero_is_entirely_about_file_conformity() -> None:
    labels = {c.control_id for c in requested(ValidationLevel.FILE_CONFORMITY)}
    assert "inventory.net_charge" in labels
    assert "torsion.profile_vs_quantum" not in labels


def test_level_zero_controls_are_available_except_geometry() -> None:
    """Tout le niveau 0 est outillé aujourd'hui, sauf la géométrie (étape 4)."""
    missing = {c.control_id for c in unavailable(ValidationLevel.FILE_CONFORMITY)}
    assert missing == {"structure.geometry"}


def test_a_requested_but_unexecuted_control_fails() -> None:
    results = audit(ValidationLevel.FILE_CONFORMITY, _produced("inventory.net_charge"))
    failed = {r.check_id for r in failures(results)}
    assert "inventory.net_charge" not in failed, "le contrôle exécuté garde son résultat"
    assert "inventory.cmap_terms" in failed
    message = next(r.message for r in results if r.check_id == "inventory.cmap_terms")
    assert "demandé mais non exécuté" in message


def test_an_unimplemented_control_fails_and_names_the_step() -> None:
    results = audit(ValidationLevel.ENERGETIC_COHERENCE, [])
    message = next(r.message for r in results if r.check_id == "torsion.profile_vs_quantum")
    assert "pas encore implémenté" in message
    assert "11D" in message


def test_audit_returns_one_result_per_requested_control() -> None:
    for level in ValidationLevel:
        results = audit(level, [])
        assert len(results) == len(requested(level))
        assert all(r.status is Status.FAIL for r in results)


def test_a_higher_level_never_hides_a_lower_level_failure() -> None:
    produced = [
        CheckResult(
            check_id="inventory.net_charge", status=Status.FAIL, message="charge non entière"
        )
    ]
    results = audit(ValidationLevel.ENERGETIC_COHERENCE, produced)
    assert any(r.check_id == "inventory.net_charge" and r.status is Status.FAIL for r in results)


def test_audit_records_the_requested_level_in_the_details() -> None:
    results = audit(ValidationLevel.NUMERIC_COHERENCE, [])
    details = next(r.details for r in results if r.check_id == "energy.single_point")
    assert details["niveau_demandé"] == "1"
    assert details["niveau_du_contrôle"] == "1"
