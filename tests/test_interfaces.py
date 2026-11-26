"""Vérifie que les interfaces d'adaptateur restent implémentables."""

from pathlib import Path

from pepforge.config import Engine, PrepConfig
from pepforge.core.peptide import Peptide
from pepforge.engines import EngineWriter
from pepforge.forcefields import ForceFieldAdapter, ParameterSet


class _DummyForceField:
    name = "dummy"
    supported_engines = frozenset({Engine.AMBER})

    def parametrize(self, peptide: Peptide) -> ParameterSet:
        return ParameterSet(forcefield=self.name, residue_names=[r.code for r in peptide.residues])


class _DummyWriter:
    engine = Engine.AMBER

    def write(self, parameters: ParameterSet, config: PrepConfig, outdir: Path) -> list[Path]:
        return []


def test_protocols_are_satisfied() -> None:
    assert isinstance(_DummyForceField(), ForceFieldAdapter)
    assert isinstance(_DummyWriter(), EngineWriter)
