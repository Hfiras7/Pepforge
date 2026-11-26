"""Tests de la configuration de préparation."""

from pathlib import Path

import pytest

from pepforge.config import Engine, PrepConfig, SimulationGoal
from pepforge.exceptions import ConfigurationError


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
