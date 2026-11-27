"""Tests de la configuration de préparation."""

from pathlib import Path

import pytest

from pepforge.config import (
    SEQUENCE_NOTATIONS,
    Engine,
    PrepConfig,
    SimulationGoal,
    ValidationLevel,
)
from pepforge.exceptions import ConfigurationError
from pepforge.io.sequence import SequenceNotation


def test_empty_config_lists_all_open_decisions() -> None:
    missing = PrepConfig().missing_decisions()
    assert "moteur de simulation" in missing
    assert "champ de force" in missing
    assert any("objectif" in m for m in missing)


def test_no_engine_or_forcefield_by_default() -> None:
    config = PrepConfig()
    assert config.engine is None
    assert config.forcefield is None
    assert config.goal is None
    with pytest.raises(ConfigurationError):
        config.require_complete()


def test_complex_requires_receptor() -> None:
    config = PrepConfig(
        sequence="Ac-RGDS-NH2",
        sequence_notation="one_letter",
        engine=Engine.AMBER,
        forcefield="ff19SB",
        water_model="OPC",
        goal=SimulationGoal.PEPTIDE_RECEPTOR_COMPLEX,
        ph=7.4,
    )
    assert config.missing_decisions() == ["récepteur (structure ou code PDB)"]
    config.receptor = "DCC"
    config.require_complete()


def test_yaml_round_trip(tmp_path: Path) -> None:
    config = PrepConfig(sequence="RGDS", engine=Engine.GROMACS, goal=SimulationGoal.FREE_PEPTIDE)
    path = tmp_path / "prep.yaml"
    config.to_yaml(path)
    assert PrepConfig.from_yaml(path) == config


def test_unknown_keys_and_values_are_rejected() -> None:
    with pytest.raises(ConfigurationError):
        PrepConfig.from_dict({"moteur": "amber"})
    with pytest.raises(ConfigurationError):
        PrepConfig.from_dict({"engine": "openmm"})


def test_solvation_has_no_default_and_round_trips(tmp_path: Path) -> None:
    config = PrepConfig()
    assert len(config.solvation.missing_decisions()) == 4
    config = PrepConfig.from_dict(
        {"solvation": {"box_shape": "truncated_octahedron", "margin_angstrom": 12.0}}
    )
    assert config.solvation.missing_decisions() == ["concentration en sel", "neutralisation"]
    path = tmp_path / "prep.yaml"
    config.to_yaml(path)
    assert PrepConfig.from_yaml(path) == config


def test_unknown_solvation_key_is_rejected() -> None:
    with pytest.raises(ConfigurationError):
        PrepConfig.from_dict({"solvation": {"padding": 10}})


def test_validation_level_defaults_to_the_free_level() -> None:
    """Le niveau 0 ne coûte rien : c'est le seul choix pris par défaut."""
    assert PrepConfig().validation_level is ValidationLevel.FILE_CONFORMITY


def test_validation_level_survives_a_round_trip(tmp_path: Path) -> None:
    config = PrepConfig(validation_level=ValidationLevel.ENERGETIC_COHERENCE)
    path = tmp_path / "prep.yaml"
    config.to_yaml(path)
    assert PrepConfig.from_yaml(path).validation_level is ValidationLevel.ENERGETIC_COHERENCE


def test_an_unknown_validation_level_is_refused() -> None:
    with pytest.raises(ConfigurationError):
        PrepConfig.from_dict({"validation_level": 9})


# --- Mode de lecture de la séquence (revue du module 1) --------------------------


def test_a_sequence_without_its_reading_mode_is_an_open_decision() -> None:
    missing = PrepConfig(sequence="ALA").missing_decisions()
    assert "mode de lecture de la séquence (code à une lettre ou à trois lettres)" in missing
    complete = PrepConfig(sequence="ALA", sequence_notation="three_letter").missing_decisions()
    assert not any("mode de lecture" in m for m in complete)


def test_no_reading_mode_is_asked_without_a_sequence() -> None:
    assert not any("mode de lecture" in m for m in PrepConfig().missing_decisions())


def test_an_unknown_reading_mode_is_refused() -> None:
    with pytest.raises(ConfigurationError, match="Mode de lecture inconnu"):
        PrepConfig(sequence="RGDS", sequence_notation="auto")
    with pytest.raises(ConfigurationError, match="Mode de lecture inconnu"):
        PrepConfig.from_dict({"sequence": "RGDS", "sequence_notation": 1})


def test_the_reading_mode_survives_a_round_trip(tmp_path: Path) -> None:
    config = PrepConfig(sequence="Ac-ALA-NHMe", sequence_notation="three_letter")
    path = tmp_path / "prep.yaml"
    config.to_yaml(path)
    assert PrepConfig.from_yaml(path).sequence_notation == "three_letter"


def test_config_and_parser_accept_the_same_reading_modes() -> None:
    """`config` ne peut pas importer `io` : les deux listes sont tenues égales ici."""
    assert tuple(n.value for n in SequenceNotation) == SEQUENCE_NOTATIONS
