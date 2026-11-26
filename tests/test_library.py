"""Tests de la bibliothèque de résidus."""

from pepforge.residues import load_default_library


def test_twenty_standard_amino_acids() -> None:
    lib = load_default_library()
    assert len(lib.residues) == 20
    assert {r.one_letter for r in lib.residues.values()} == set("ACDEFGHIKLMNPQRSTVWY")


def test_every_data_file_declares_its_source() -> None:
    lib = load_default_library()
    assert lib.sources
    assert all(source.strip() for source in lib.sources.values())
