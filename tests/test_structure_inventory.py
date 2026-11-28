"""Tests de l'inventaire d'une structure et de la correspondance `Structure` → `Peptide`."""

from __future__ import annotations

from pathlib import Path

import pytest

from pepforge.core.structure import (
    Atom,
    Structure,
)
from pepforge.io import inventory, read_structure
from pepforge.report.journal import Journal
from pepforge.residues.structure_components import load_structure_components

DATA = Path(__file__).parent / "data" / "structures"


def _read(name: str, model: int | None = None) -> Structure:
    return read_structure(DATA / name, journal=Journal(), model=model)


def _residue(
    chain: str, number: int, name: str, atoms: tuple[str, ...], icode: str = ""
) -> list[Atom]:
    return [
        Atom(
            i,
            atom,
            atom[0],
            name,
            number,
            chain,
            float(number),
            float(i),
            0.0,
            insertion_code=icode,
        )
        for i, atom in enumerate(atoms, start=1)
    ]


BACKBONE = ("N", "CA", "C", "O")


def _chain(*residues: tuple[int, str]) -> tuple[Atom, ...]:
    atoms: list[Atom] = []
    for number, name in residues:
        atoms += _residue("A", number, name, BACKBONE)
    return tuple(atoms)


# --- Données de référence ---------------------------------------------------------------


def test_structure_components_are_sourced_and_consistent() -> None:
    components = load_structure_components()
    assert "Chemical Component Dictionary" in components.source
    assert components.water == {"HOH", "DOD"}
    assert components.n_terminal_caps == {"ACE": "acetyl"}
    assert components.c_terminal_caps == {"NH2": "amide", "NME": "n_methylamide"}
    assert components.free_c_terminus == "free_acid"
    assert components.free_n_terminus == "free_amine"


# --- Inventaire -------------------------------------------------------------------------


def test_inventory_of_the_capped_peptide() -> None:
    inv = inventory(_read("capped_peptide.pdb"))
    assert inv.format == "pdb"
    assert inv.atoms == 32
    assert inv.hydrogens == 0
    (chain,) = inv.chains
    assert chain.chain_id == "A"
    assert chain.standard_residues == 5
    assert chain.waters == 1
    assert chain.ions == ("A:ZN201",)
    assert chain.unrecognized == {"ACE": 1, "NH2": 1}
    assert chain.insertion_codes == ("A:GLY3A",)
    assert chain.numbering_gaps == ()
    assert chain.first_residue == "A:ALA1" and chain.last_residue == "A:CYS4"
    (site,) = inv.alternate_locations
    assert site.residue.label() == "A:ALA1"
    assert site.alt_ids == ("A", "B")
    assert site.atoms == 2
    assert site.residue_names == ()
    assert inv.connections == {"disulfide": 1}


def test_inventory_lists_every_decision_to_take() -> None:
    decisions = " | ".join(inventory(_read("capped_peptide.pdb")).pending_decisions)
    for expected in (
        "conformations alternatives",
        "Composants non reconnus",
        "ACE×1",
        "1 ion(s)",
        "1 molécule(s) d'eau",
        "Aucun hydrogène",
        "Liaisons déclarées",
    ):
        assert expected in decisions


def test_inventory_of_an_ensemble_names_the_model_read() -> None:
    inv = inventory(_read("nmr_two_models.pdb", model=2))
    assert inv.model == 2 and inv.model_count == 2
    assert inv.hydrogens == 1
    assert any("seul le modèle 2 est lu" in d for d in inv.pending_decisions)


def test_numbering_gap_is_reported_not_repaired() -> None:
    structure = Structure(atoms=_chain((1, "ALA"), (2, "GLY"), (5, "SER")))
    inv = inventory(structure)
    assert inv.chains[0].numbering_gaps == (("A:GLY2", "A:SER5"),)
    assert any("saut(s) de numérotation" in d for d in inv.pending_decisions)
    assert len(structure.atoms) == 12  # rien n'a été ajouté ni retiré


def test_modified_amino_acid_counts_for_numbering_but_is_unrecognized() -> None:
    # MSE (sélénométhionine) porte N, CA et C : il compte dans la chaîne, sans être standard.
    structure = Structure(atoms=_chain((1, "ALA"), (2, "MSE"), (3, "GLY")))
    chain = inventory(structure).chains[0]
    assert chain.amino_acid_like_residues == 3
    assert chain.standard_residues == 2
    assert chain.unrecognized == {"MSE": 1}
    assert chain.numbering_gaps == ()


def test_microheterogeneity_is_reported() -> None:
    first = [
        Atom(1, "CA", "C", "SER", 7, "A", 0.0, 0.0, 0.0, alt_loc="A", occupancy=0.5),
        Atom(2, "CA", "C", "THR", 7, "A", 0.1, 0.0, 0.0, alt_loc="B", occupancy=0.5),
    ]
    inv = inventory(Structure(atoms=tuple(first)))
    (site,) = inv.alternate_locations
    assert site.residue_names == ("SER", "THR")
    assert any(
        "change selon la conformation" in d or "nom change" in d for d in inv.pending_decisions
    )


def test_zero_and_missing_occupancies_are_counted() -> None:
    atoms = (
        Atom(1, "CA", "C", "ALA", 1, "A", 0.0, 0.0, 0.0, occupancy=0.0),
        Atom(2, "CB", "C", "ALA", 1, "A", 1.0, 0.0, 0.0),
    )
    inv = inventory(Structure(atoms=atoms))
    assert inv.atoms_with_zero_occupancy == 1
    assert inv.atoms_without_occupancy == 1


def test_inventory_is_serializable() -> None:
    import yaml

    data = inventory(_read("capped_peptide.cif")).to_dict()
    assert yaml.safe_load(yaml.safe_dump(data, allow_unicode=True))["atoms"] == 32


def test_peptide_bond_criterion_is_sourced_and_can_be_imposed() -> None:
    from pepforge.exceptions import ReferenceDataError
    from pepforge.residues.structure_components import load_peptide_bond_criterion

    criterion = load_peptide_bond_criterion()
    assert criterion.max_c_n_distance == 1.8
    assert "PPBuilder" in criterion.source
    assert criterion.reference_length == 1.329
    assert "Engh" in criterion.reference_source
    imposed = load_peptide_bond_criterion(2.2)
    assert imposed.max_c_n_distance == 2.2
    assert "imposée" in imposed.source and "1.8" in imposed.source
    with pytest.raises(ReferenceDataError, match="inférieur ou égal"):
        load_peptide_bond_criterion(1.3)
