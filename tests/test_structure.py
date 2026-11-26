"""Tests du modèle de structure 3D."""

from __future__ import annotations

import dataclasses

import pytest

from pepforge.core import Atom, Structure


def _atom(serial: int, element: str = "C", chain: str = "A", alt: str | None = None) -> Atom:
    return Atom(serial, "X", element, "ALA", 1, chain, 0.0, 0.0, 0.0, alt_loc=alt)


def test_structure_summary() -> None:
    structure = Structure(atoms=(_atom(1), _atom(2, "H", "B"), _atom(3, alt="A")))
    assert len(structure) == 3
    assert structure.chain_ids == ["A", "B"]
    assert structure.has_hydrogens()
    assert structure.has_alternate_locations()


def test_structure_is_immutable() -> None:
    structure = Structure(atoms=(_atom(1),))
    with pytest.raises(dataclasses.FrozenInstanceError):
        structure.source = "modifié"  # type: ignore[misc]
