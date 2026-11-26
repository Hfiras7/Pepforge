"""Tests de contrat réutilisables.

Chaque implémentation réelle d'une interface (adaptateur de champ de force,
générateur moteur, constructeur 3D, étape) doit passer la fonction de contrat
correspondante, en plus de ses propres tests. Exemple, pour un futur
adaptateur :

    def test_ff19sb_contract() -> None:
        assert_forcefield_contract(FF19SBAdapter(), [parse_sequence("Ac-RGDS-NH2")])
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from pepforge.build import StructureBuilder
from pepforge.config import PrepConfig
from pepforge.core import Peptide, Structure
from pepforge.engines import EngineWriter
from pepforge.exceptions import MissingParameterError
from pepforge.forcefields import ForceFieldAdapter, ParameterSet
from pepforge.pipeline import Stage, StageContext
from pepforge.report import Journal


def assert_forcefield_contract(adapter: ForceFieldAdapter, peptides: Sequence[Peptide]) -> None:
    assert isinstance(adapter, ForceFieldAdapter)
    assert adapter.supported_engines, "un adaptateur déclare au moins un moteur"
    for peptide in peptides:
        params = adapter.parametrize(peptide)
        assert isinstance(params, ParameterSet)
        assert params.forcefield == adapter.name
        assert params.missing == [], "paramètres manquants : lever MissingParameterError"
        assert params.residue_names, "aucun résidu paramétré"
        assert params.sources, "la provenance des paramètres est obligatoire"


def assert_engine_writer_contract(
    writer: EngineWriter, structure: Structure, tmp_path: Path
) -> None:
    assert isinstance(writer, EngineWriter)
    incomplete = ParameterSet(forcefield="x", residue_names=["ALA"], missing=["CT-N-XX"])
    with pytest.raises(MissingParameterError):
        writer.write_topology(structure, incomplete, tmp_path)


def assert_builder_contract(builder: StructureBuilder, peptide: Peptide) -> None:
    assert isinstance(builder, StructureBuilder)
    journal = Journal()
    structure = builder.build(peptide, journal)
    assert isinstance(structure, Structure)
    assert structure.peptide is not None
    assert structure.peptide.n_terminus == peptide.n_terminus, "extrémité N perdue"
    assert structure.peptide.c_terminus == peptide.c_terminus, "extrémité C perdue"
    assert journal.entries, "la construction doit être journalisée"


def assert_stage_contract(stage: Stage[object, object], data: object, tmp_path: Path) -> None:
    assert isinstance(stage, Stage)
    assert stage.name, "une étape porte un nom stable (utilisé dans le rapport)"
    context = StageContext(config=PrepConfig(), workdir=tmp_path)
    stage.run(data, context)
