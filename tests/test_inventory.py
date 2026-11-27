"""Tests de la comparaison d'inventaires, cœur de l'absence d'erreur silencieuse.

Chaque test correspond à une ligne du catalogue de `docs/ARCHITECTURE.md`,
section 17.2 : un piège qui, sans ce contrôle, passerait inaperçu.
"""

from __future__ import annotations

from pepforge.checks import CheckResult, Inventory, Status, compare, expected_from_peptide
from pepforge.checks.base import failures
from pepforge.io import parse_sequence


def _complete(**overrides: object) -> Inventory:
    """Inventaire entièrement renseigné, donc intégralement vérifiable."""
    base: dict[str, object] = {
        "residue_sequence": ("ARG", "GLY", "ASP", "SER"),
        "n_terminus": "acetyl",
        "c_terminus": "amide",
        "atom_counts": {"1:ARG": 24, "2:GLY": 7, "3:ASP": 12, "4:SER": 11},
        "backbone_bonds": 3,
        "disulfides": (),
        "net_charge": 0.0,
        "cmap_terms": 4,
        "source": "test",
    }
    base.update(overrides)
    return Inventory(**base)  # type: ignore[arg-type]


def _result(results: list[CheckResult], check_id: str) -> CheckResult:
    for result in results:
        if result.check_id == check_id:
            return result
    raise AssertionError(f"contrôle absent du résultat : {check_id}")


def _status(results: list[CheckResult], check_id: str) -> Status:
    return _result(results, check_id).status


def test_identical_inventories_all_pass() -> None:
    results = compare(_complete(), _complete())
    assert failures(results) == []
    assert len(results) == 8, "un résultat par champ, même quand tout va bien"


def test_dropped_residue_is_caught() -> None:
    """Piège 17.2 : un résidu ignoré faute de topologie, le fichier sort amputé."""
    produced = _complete(
        residue_sequence=("ARG", "GLY", "SER"),
        atom_counts={"1:ARG": 24, "2:GLY": 7, "4:SER": 11},
        backbone_bonds=2,
    )
    results = compare(_complete(), produced)
    assert _status(results, "inventory.residue_sequence") is Status.FAIL
    assert _status(results, "inventory.atom_counts") is Status.FAIL
    assert _status(results, "inventory.backbone_bonds") is Status.FAIL


def test_substituted_residue_reports_its_position() -> None:
    produced = _complete(residue_sequence=("ARG", "GLY", "ASN", "SER"))
    results = compare(_complete(), produced)
    message = _result(results, "inventory.residue_sequence").message
    assert _status(results, "inventory.residue_sequence") is Status.FAIL
    assert "position 3" in message
    assert "ASP" in message and "ASN" in message


def test_missing_atoms_in_one_residue() -> None:
    produced = _complete(atom_counts={"1:ARG": 24, "2:GLY": 7, "3:ASP": 10, "4:SER": 11})
    results = compare(_complete(), produced)
    assert _status(results, "inventory.atom_counts") is Status.FAIL


def test_wrong_terminus_is_caught() -> None:
    """Piège 17.2 : patch non appliqué, donc extrémité fausse."""
    results = compare(_complete(), _complete(c_terminus="free_acid"))
    assert _status(results, "inventory.c_terminus") is Status.FAIL


def test_chain_break_is_caught() -> None:
    """Piège 17.2 : résidus manquants, donc rupture du squelette."""
    results = compare(_complete(), _complete(backbone_bonds=2))
    assert _status(results, "inventory.backbone_bonds") is Status.FAIL


def test_missing_disulfide_is_caught() -> None:
    expected = _complete(disulfides=((2, 7),))
    results = compare(expected, _complete(disulfides=()))
    assert _status(results, "inventory.disulfides") is Status.FAIL


def test_omitted_cmap_term_is_caught() -> None:
    """Piège 17.2 : un terme CMAP omis parce qu'une référence ne s'est pas résolue."""
    results = compare(_complete(), _complete(cmap_terms=3))
    assert _status(results, "inventory.cmap_terms") is Status.FAIL


def test_non_integral_charge_is_caught() -> None:
    """Piège 17.2 : redistribution de jonction approximative."""
    results = compare(_complete(), _complete(net_charge=0.02))
    assert _status(results, "inventory.net_charge") is Status.FAIL
    assert "tolérance" in _result(results, "inventory.net_charge").message


def test_rounding_noise_on_charge_passes() -> None:
    results = compare(_complete(), _complete(net_charge=1e-9))
    assert _status(results, "inventory.net_charge") is Status.PASS


def test_unverifiable_field_fails_and_never_passes() -> None:
    """La décision centrale : « je n'ai pas pu contrôler » n'est pas « tout va bien »."""
    results = compare(_complete(), _complete(cmap_terms=None))
    assert _status(results, "inventory.cmap_terms") is Status.FAIL
    message = _result(results, "inventory.cmap_terms").message
    assert "non vérifiable" in message
    assert "fichier produit" in message


def test_unverifiable_on_the_expected_side_also_fails() -> None:
    results = compare(_complete(net_charge=None), _complete())
    assert _status(results, "inventory.net_charge") is Status.FAIL
    message = _result(results, "inventory.net_charge").message
    assert "non vérifiable" in message
    assert "attendu" in message


def test_comparison_does_not_stop_at_the_first_difference() -> None:
    produced = _complete(c_terminus="free_acid", cmap_terms=3, net_charge=0.5)
    failed = {r.check_id for r in failures(compare(_complete(), produced))}
    assert failed == {
        "inventory.c_terminus",
        "inventory.cmap_terms",
        "inventory.net_charge",
    }


def test_results_carry_both_sources_for_the_report() -> None:
    results = compare(_complete(source="attendu"), _complete(source="psf relu"))
    details = results[0].details
    assert details["source_attendu"] == "attendu"
    assert details["source_produit"] == "psf relu"


def test_expected_from_peptide_fills_what_is_deducible() -> None:
    inventory = expected_from_peptide(parse_sequence("Ac-RGDS-NH2"))
    assert inventory.residue_sequence == ("ARG", "GLY", "ASP", "SER")
    assert inventory.n_terminus == "acetyl"
    assert inventory.c_terminus == "amide"
    assert inventory.backbone_bonds == 3
    # Ces champs dépendent de la bibliothèque de résidus et du champ de force,
    # qui ne sont pas encore écrits : ils restent inconnus, donc non vérifiables.
    assert inventory.atom_counts is None
    assert inventory.net_charge is None
    assert inventory.cmap_terms is None


def test_incomplete_expected_inventory_blocks_rather_than_passes() -> None:
    """Un inventaire attendu encore incomplet ne laisse rien passer en silence."""
    results = compare(expected_from_peptide(parse_sequence("Ac-RGDS-NH2")), _complete())
    failed = {r.check_id for r in failures(results)}
    assert failed == {
        "inventory.atom_counts",
        "inventory.net_charge",
        "inventory.cmap_terms",
    }
