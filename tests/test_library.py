"""Tests de la bibliothèque de résidus."""

import pytest

from pepforge.exceptions import UnknownTerminusError
from pepforge.residues import load_default_library


def test_twenty_standard_amino_acids() -> None:
    lib = load_default_library()
    assert len(lib.residues) == 20
    assert {r.one_letter for r in lib.residues.values()} == set("ACDEFGHIKLMNPQRSTVWY")


def test_every_data_file_declares_its_source() -> None:
    lib = load_default_library()
    assert lib.sources
    assert all(source.strip() for source in lib.sources.values())


def test_each_terminus_is_written_with_its_first_alias() -> None:
    lib = load_default_library()
    assert lib.n_terminus_symbol("acetyl") == "Ac"
    assert lib.n_terminus_symbol("free_amine") == "H"
    assert lib.c_terminus_symbol("amide") == "NH2"
    assert lib.c_terminus_symbol("free_acid") == "OH"
    assert lib.c_terminus_symbol("n_methylamide") == "NHMe"


def test_every_written_terminus_reads_back_as_itself() -> None:
    lib = load_default_library()
    assert set(lib.n_symbols) == set(lib.n_termini)
    assert set(lib.c_symbols) == set(lib.c_termini)
    for key in lib.n_termini:
        assert lib.n_terminus(lib.n_terminus_symbol(key)).key == key
    for key in lib.c_termini:
        assert lib.c_terminus(lib.c_terminus_symbol(key)).key == key


def test_an_unknown_terminus_has_no_written_form() -> None:
    lib = load_default_library()
    with pytest.raises(UnknownTerminusError):
        lib.n_terminus_symbol("formyl")
    with pytest.raises(UnknownTerminusError):
        lib.c_terminus_symbol("thioester")
