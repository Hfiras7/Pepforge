"""Vérifie que les interfaces restent implémentables et que les contrats sont testables.

Les implémentations factices ci-dessous sont minimales : elles montrent la
forme attendue d'un adaptateur, d'un générateur moteur et d'un constructeur.
"""

from __future__ import annotations

from pathlib import Path

from contracts import (
    assert_builder_contract,
    assert_engine_writer_contract,
    assert_forcefield_contract,
    assert_stage_contract,
)
from pepforge.config import Engine, SolvationSpec
from pepforge.core import Atom, Peptide, Structure
from pepforge.engines import EngineFiles
from pepforge.exceptions import MissingParameterError
from pepforge.forcefields import ParameterSet
from pepforge.io import SequenceNotation, parse_sequence
from pepforge.pipeline import SequenceInputStage
from pepforge.report import Journal


class _DummyForceField:
    name = "dummy"
    supported_engines = frozenset({Engine.AMBER})
    compatible_water_models = frozenset({"dummy-water"})

    def parametrize(self, peptide: Peptide) -> ParameterSet:
        return ParameterSet(
            forcefield=self.name,
            residue_names=[r.code for r in peptide.residues],
            sources={"dummy": "jeu factice pour les tests"},
        )


class _DummyWriter:
    engine = Engine.AMBER

    def write_topology(
        self, structure: Structure, parameters: ParameterSet, outdir: Path
    ) -> EngineFiles:
        if parameters.missing:
            raise MissingParameterError(", ".join(parameters.missing))
        return EngineFiles(engine=self.engine, files={}, atom_count=len(structure))

    def solvate(
        self, files: EngineFiles, spec: SolvationSpec, outdir: Path, journal: Journal
    ) -> EngineFiles:
        journal.record("solvate", "test", box_shape=spec.box_shape)
        return EngineFiles(engine=self.engine, files=files.files, solvated=True)


class _DummyBuilder:
    name = "dummy-builder"

    def build(self, peptide: Peptide, journal: Journal) -> Structure:
        atoms = tuple(
            Atom(i + 1, "CA", "C", r.code, i + 1, "A", float(i), 0.0, 0.0)
            for i, r in enumerate(peptide.residues)
        )
        journal.record("build", "test", builder=self.name)
        return Structure(atoms=atoms, peptide=peptide, source=self.name)


def test_forcefield_contract() -> None:
    assert_forcefield_contract(
        _DummyForceField(), [parse_sequence("Ac-RGDS-NH2", notation="one_letter")]
    )


def test_engine_writer_contract(tmp_path: Path) -> None:
    structure = _DummyBuilder().build(
        parse_sequence("Ac-RGDS-NH2", notation="one_letter"), Journal()
    )
    assert_engine_writer_contract(_DummyWriter(), structure, tmp_path)


def test_builder_contract() -> None:
    assert_builder_contract(_DummyBuilder(), parse_sequence("Ac-RGDS-NH2", notation="one_letter"))


def test_sequence_stage_contract(tmp_path: Path) -> None:
    assert_stage_contract(SequenceInputStage(SequenceNotation.ONE_LETTER), "Ac-RGDS-NH2", tmp_path)
