"""Tests de l'étape de protonation, de ses contrôles et des charges de référence."""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pytest
import yaml

from pepforge import exceptions
from pepforge.checks import PROTONATION_CHECKS, Status, run_checks
from pepforge.config import PrepConfig
from pepforge.core import GroupAssignment, IonizationForm, Link, LinkKind, ProtonationResult
from pepforge.exceptions import ConfigurationError, QualityCheckError
from pepforge.io import SequenceNotation
from pepforge.io.sequence import parse_sequence
from pepforge.pipeline import ProtonationStage, SequenceInputStage, StageContext, run_stage
from pepforge.residues import assign_protonation

REFERENCE = Path(__file__).parent / "reference" / "charges.yaml"


def _context(tmp_path: Path, **config: Any) -> StageContext:
    return StageContext(config=PrepConfig(**config), workdir=tmp_path)


def _results(result: ProtonationResult) -> dict[str, Any]:
    return {r.check_id: r for r in run_checks(PROTONATION_CHECKS, result)}


# --- Étape --------------------------------------------------------------------


def test_missing_ph_raises_configuration_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="pH"):
        ProtonationStage().run(
            parse_sequence("H-GGG-OH", notation="one_letter"), _context(tmp_path)
        )


def test_sequence_then_protonation_through_run_stage(tmp_path: Path) -> None:
    """Intégration : `SequenceInputStage` puis `ProtonationStage`, contrôles compris."""
    context = _context(tmp_path, ph=7.4)
    peptide = run_stage(SequenceInputStage(SequenceNotation.ONE_LETTER), "Ac-KDEY-NH2", context)
    result = run_stage(ProtonationStage(), peptide, context, PROTONATION_CHECKS)
    assert result.net_charge == -1
    statuses = {r.check_id: r.status for r in context.check_results["normalize.protonation"]}
    assert set(statuses.values()) == {Status.PASS}


def test_each_choice_is_journaled_with_its_source(tmp_path: Path) -> None:
    context = _context(tmp_path, ph=7.4)
    ProtonationStage().run(parse_sequence("H-KG-OH", notation="one_letter"), context)
    entries = [e for e in context.journal.entries if e.action == "assign_protonation"]
    assert [e.parameters["site"] for e in entries] == ["N", "1", "C"]
    for entry in entries:
        assert "Grimsley" in entry.parameters["pka_source"]
        assert entry.parameters["ph"] == 7.4
        assert "pka" in entry.parameters
    (charge,) = [e for e in context.journal.entries if e.action == "net_charge"]
    assert charge.parameters["net_charge"] == 1


def test_configured_states_are_imposed_and_journaled(tmp_path: Path) -> None:
    context = _context(tmp_path, ph=7.4, protonation={"2": "his_epsilon"})
    result = ProtonationStage().run(parse_sequence("Ac-GHG-NH2", notation="one_letter"), context)
    assert result.peptide.residue_at(2).protonation == "his_epsilon"
    assert result.assignments[0].imposed
    actions = [e.action for e in context.journal.entries]
    assert actions.index("impose_protonation") < actions.index("assign_protonation")


def test_configured_terminus_state(tmp_path: Path) -> None:
    context = _context(tmp_path, ph=7.4, protonation={"N": "nterm_deprotonated"})
    result = ProtonationStage().run(parse_sequence("H-GGG-OH", notation="one_letter"), context)
    assert result.peptide.n_terminus.protonation == "nterm_deprotonated"
    assert result.net_charge == -1


@pytest.mark.parametrize("site", ["0", "4", "X", "-1", "n"])
def test_invalid_configured_site_is_refused(tmp_path: Path, site: str) -> None:
    context = _context(tmp_path, ph=7.4, protonation={site: "his_epsilon"})
    with pytest.raises(ConfigurationError, match="site"):
        ProtonationStage().run(parse_sequence("Ac-GHG-NH2", notation="one_letter"), context)


def test_protonation_survives_a_config_round_trip(tmp_path: Path) -> None:
    config = PrepConfig(ph=7.4, protonation={"2": "his_epsilon", "N": "nterm_protonated"})
    path = tmp_path / "config.yaml"
    config.to_yaml(path)
    assert PrepConfig.from_yaml(path).protonation == config.protonation


def test_integer_site_keys_from_yaml_are_normalized() -> None:
    config = PrepConfig.from_dict({"ph": 7.4, "protonation": {2: "his_epsilon"}})
    assert config.protonation == {"2": "his_epsilon"}


def test_protonation_must_be_a_mapping() -> None:
    with pytest.raises(ConfigurationError):
        PrepConfig.from_dict({"protonation": ["his_epsilon"]})


def test_near_pka_is_a_warning_not_a_failure(tmp_path: Path) -> None:
    context = _context(tmp_path, ph=7.4)
    run_stage(
        ProtonationStage(),
        parse_sequence("H-GGG-OH", notation="one_letter"),
        context,
        PROTONATION_CHECKS,
    )
    (entry,) = [e for e in context.journal.entries if e.action == "stage_completed"]
    assert entry.parameters["warnings"] == ["protonation.near_pka"]


# --- Contrôles ----------------------------------------------------------------


def _assignment(result: ProtonationResult, index: int = 0) -> GroupAssignment:
    return result.assignments[index]


def test_checks_pass_on_a_clean_result() -> None:
    result = assign_protonation(parse_sequence("Ac-KDEY-NH2", notation="one_letter"), 7.4)
    assert {r.status for r in _results(result).values()} == {Status.PASS}


def test_state_missing_on_the_peptide_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    broken = dataclasses.replace(result, peptide=result.peptide.with_protonation(2, None))
    assert _results(broken)["protonation.states_recorded"].status is Status.FAIL


def test_stray_state_on_the_peptide_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    broken = dataclasses.replace(result, peptide=result.peptide.with_protonation(1, "his_delta"))
    assert _results(broken)["protonation.states_recorded"].status is Status.FAIL


def test_wrong_net_charge_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    assert _results(dataclasses.replace(result, net_charge=0))["protonation.net_charge"].status is (
        Status.FAIL
    )


def test_non_integer_net_charge_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    broken = dataclasses.replace(result, net_charge=1.0)
    assert _results(broken)["protonation.net_charge"].status is Status.FAIL


def test_unsourced_pka_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    bad = dataclasses.replace(_assignment(result), pka_source=" ")
    broken = dataclasses.replace(result, assignments=(bad,))
    assert _results(broken)["protonation.sources"].status is Status.FAIL


def test_deduced_minority_state_fails() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), 7.4)
    bad = dataclasses.replace(_assignment(result), form=IonizationForm.BASE)
    broken = dataclasses.replace(result, assignments=(bad,))
    assert _results(broken)["protonation.deduced_states_follow_ph"].status is Status.FAIL


def test_imposed_minority_state_warns() -> None:
    peptide = parse_sequence("Ac-GDG-NH2", notation="one_letter").with_protonation(
        2, "asp_protonated"
    )
    results = _results(assign_protonation(peptide, 7.4))
    assert results["protonation.imposed_minority"].status is Status.WARN
    assert results["protonation.deduced_states_follow_ph"].status is Status.PASS


def test_imposed_state_at_pka_warns() -> None:
    peptide = parse_sequence("Ac-GKG-NH2", notation="one_letter").with_protonation(
        2, "lys_protonated"
    )
    results = _results(assign_protonation(peptide, 10.4))
    assert results["protonation.imposed_minority"].status is Status.WARN
    assert results["protonation.near_pka"].status is Status.WARN


@pytest.mark.parametrize(("ph", "status"), [(9.41, Status.WARN), (9.4, Status.PASS)])
def test_near_pka_window_boundary(ph: float, status: Status) -> None:
    """Écart strictement inférieur à 1 unité : alerte ; écart de 1 : pas d'alerte."""
    results = _results(assign_protonation(parse_sequence("Ac-GKG-NH2", notation="one_letter"), ph))
    assert results["protonation.near_pka"].status is status


def test_stage_fails_when_a_check_fails(tmp_path: Path) -> None:
    class _Broken(ProtonationStage):
        def run(self, data: Any, context: StageContext) -> ProtonationResult:
            result = super().run(data, context)
            return dataclasses.replace(result, net_charge=result.net_charge + 1)

    with pytest.raises(QualityCheckError, match=r"protonation\.net_charge"):
        run_stage(
            _Broken(),
            parse_sequence("Ac-GKG-NH2", notation="one_letter"),
            _context(tmp_path, ph=7.4),
            PROTONATION_CHECKS,
        )


# --- Validation : charges de référence calculées à la main --------------------


def _reference() -> dict[str, Any]:
    data = yaml.safe_load(REFERENCE.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _ids(entries: list[dict[str, Any]]) -> list[str]:
    return [f"{e['sequence']}@{e['ph']}" for e in entries]


@pytest.mark.parametrize("entry", _reference()["peptides"], ids=_ids(_reference()["peptides"]))
def test_reference_charges(entry: dict[str, Any], tmp_path: Path) -> None:
    peptide = parse_sequence(entry["sequence"], notation="one_letter")
    if "disulfide" in entry:
        peptide = peptide.with_link(Link(LinkKind.DISULFIDE, *entry["disulfide"]))
    context = _context(tmp_path, ph=entry["ph"], protonation=entry.get("imposed", {}))
    result = run_stage(ProtonationStage(), peptide, context, PROTONATION_CHECKS)
    assert result.net_charge == entry["net_charge"]
    tolerance = _reference()["tolerance_theoretical"]
    assert result.theoretical_net_charge == pytest.approx(entry["theoretical"], abs=tolerance)
    near = _results(result)["protonation.near_pka"]
    flagged = [a.site for a in result.assignments if abs(result.ph - a.pka) < 1.0]
    assert flagged == entry["near_pka"]
    assert (near.status is Status.WARN) == bool(entry["near_pka"])


@pytest.mark.parametrize("entry", _reference()["refused"], ids=_ids(_reference()["refused"]))
def test_reference_refusals(entry: dict[str, Any], tmp_path: Path) -> None:
    error = getattr(exceptions, entry["error"])
    with pytest.raises(error):
        ProtonationStage().run(
            parse_sequence(entry["sequence"], notation="one_letter"),
            _context(tmp_path, ph=entry["ph"]),
        )
