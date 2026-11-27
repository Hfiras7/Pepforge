"""Tests de l'interface de calcul quantique.

Aucun logiciel quantique n'est installé ici : ces tests portent sur ce qui doit
être vrai sans eux, c'est-à-dire la spécification du niveau de théorie et la
provenance d'un profil. Une implémentation réelle d'ORCA ou de Psi4 passera en
plus `tests/contracts.py::assert_quantum_tool_contract`.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from contracts import assert_quantum_tool_contract
from pepforge.exceptions import ConfigurationError
from pepforge.tools import (
    QuantumSoftware,
    QuantumSpec,
    QuantumTool,
    ScanPoint,
    TorsionProfile,
)


def _spec() -> QuantumSpec:
    return QuantumSpec(
        software=QuantumSoftware.ORCA,
        method="MP2",
        basis_set="6-31G(d)",
        relaxed_scan=True,
        scan_step_degrees=10.0,
    )


class _FakeOrca:
    """Faux logiciel, uniquement pour vérifier que le contrat est vérifiable."""

    software = QuantumSoftware.ORCA

    def version(self) -> str:
        return "6.0.0 (simulée)"

    def optimize(self, geometry: Path, spec: QuantumSpec, outdir: Path) -> Path:
        return geometry

    def scan_dihedral(
        self,
        geometry: Path,
        dihedral: Sequence[int],
        spec: QuantumSpec,
        outdir: Path,
    ) -> TorsionProfile:
        points = tuple(
            ScanPoint(angle_degrees=float(i * 10), energy_hartree=-100.0 + i * 1e-4)
            for i in range(36)
        )
        return TorsionProfile(
            fragment="aza-beta3-ala-capé",
            dihedral=(dihedral[0], dihedral[1], dihedral[2], dihedral[3]),
            points=points,
            spec=spec,
            software_version=self.version(),
            computed_on="2026-11-27",
            hardware="poste de bureau, 8 cœurs",
        )

    def charges(self, geometry: Path, spec: QuantumSpec, outdir: Path) -> dict[str, float]:
        return {}


def test_an_empty_spec_lists_every_open_decision() -> None:
    """Le niveau de théorie ne se devine pas : il se déclare."""
    missing = QuantumSpec().missing_decisions()
    assert "logiciel de calcul quantique" in missing
    assert "méthode" in missing
    assert "base" in missing


def test_an_incomplete_spec_refuses_to_produce_a_label() -> None:
    with pytest.raises(ConfigurationError):
        QuantumSpec(software=QuantumSoftware.PSI4).label()


def test_a_complete_spec_labels_the_level_of_theory() -> None:
    assert _spec().label() == "MP2/6-31G(d)"


def test_relative_energies_start_at_zero_and_are_in_kcal() -> None:
    points = (
        ScanPoint(angle_degrees=0.0, energy_hartree=-100.0),
        ScanPoint(angle_degrees=10.0, energy_hartree=-100.0 + 1.0 / 627.509474063),
    )
    profile = TorsionProfile(
        fragment="f", dihedral=(1, 2, 3, 4), points=points, spec=_spec(), software_version="x"
    )
    relative = profile.relative_energies_kcal()
    assert relative[0] == 0.0
    assert relative[1] == pytest.approx(1.0)


def test_an_empty_profile_has_no_relative_energies() -> None:
    profile = TorsionProfile(
        fragment="f", dihedral=(1, 2, 3, 4), points=(), spec=_spec(), software_version="x"
    )
    assert profile.relative_energies_kcal() == ()


def test_a_fake_tool_satisfies_the_protocol_and_its_contract(tmp_path: Path) -> None:
    geometry = tmp_path / "fragment.xyz"
    geometry.write_text("", encoding="utf-8")
    tool: QuantumTool = _FakeOrca()
    assert_quantum_tool_contract(tool, _spec(), geometry)
