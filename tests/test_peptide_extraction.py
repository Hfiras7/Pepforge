"""Tests de l'extraction d'un peptide depuis une structure, et de son rapport.

Deux familles : des chaînes construites à la main, où chaque catégorie est
provoquée et sa géométrie connue, et trois structures réelles (4ZHL, 6WQA,
7DDO, `tests/data/structures/reference/`) qui ont révélé les défauts corrigés.
"""

from __future__ import annotations

import gzip
import hashlib
from functools import cache
from pathlib import Path

import pytest
import yaml

from pepforge.core.peptide import C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED, LinkKind
from pepforge.core.structure import (
    Atom,
    Connection,
    ConnectionKind,
    MissingAtom,
    ResidueId,
    Structure,
)
from pepforge.exceptions import ExtractionRefusedError, StructureError
from pepforge.io import (
    Decision,
    ExtractionDecisions,
    ExtractionReport,
    analyze_chain,
    extract_peptides,
    read_structure,
    structure_to_peptide,
)
from pepforge.io.peptide_extraction import (
    STRUCTURE_CHAIN,
    STRUCTURE_SEGMENT,
    ComponentCategory,
    JunctionCategory,
)
from pepforge.report.journal import Journal
from pepforge.residues.structure_components import load_peptide_bond_criterion

DATA = Path(__file__).parent / "data" / "structures"
REFERENCE = DATA / "reference"
BOND = 1.33
"""Distance C–N des chaînes construites à la main (proche de 1,329 Å, Engh et Huber 1991)."""

BACKBONE = ("N", "CA", "C", "O")


# --- Construction de chaînes à géométrie connue -------------------------------------------


def _build(
    *specs: tuple[int, str] | tuple[int, str, dict[str, object]],
    chain: str = "A",
) -> tuple[Atom, ...]:
    """Résidus alignés sur x ; option `gap` = distance C–N avec le précédent.

    Options d'un résidu : `icode`, `gap` (Å), `atoms` (noms d'atomes), `extra`
    (atomes ajoutés, placés à l'écart du squelette).
    """
    atoms: list[Atom] = []
    x = 0.0
    serial = 0
    for k, spec in enumerate(specs):
        number, name = spec[0], spec[1]
        options: dict[str, object] = spec[2] if len(spec) == 3 else {}  # type: ignore[misc]
        if k:
            x += float(options.get("gap", BOND))  # type: ignore[arg-type]
        names = tuple(options.get("atoms", BACKBONE))  # type: ignore[call-overload]
        names += tuple(options.get("extra", ()))  # type: ignore[call-overload]
        positions = {"N": (x, 0.0), "CA": (x + 1.0, 0.5), "C": (x + 2.0, 0.0), "O": (x + 2.0, 1.2)}
        for offset, atom in enumerate(names):
            serial += 1
            px, py = positions.get(atom, (x + 1.0, 2.0 + offset))
            atoms.append(
                Atom(
                    serial,
                    atom,
                    "H" if atom.startswith("H") else atom[0],
                    name,
                    number,
                    chain,
                    px,
                    py,
                    0.0,
                    insertion_code=str(options.get("icode", "")),
                )
            )
        x += 2.0
    return tuple(atoms)


def _rid(number: int, name: str, chain: str = "A", icode: str = "") -> ResidueId:
    return ResidueId(chain, number, icode, name)


def _refused(structure: Structure, chain: str = "A", **kwargs: object) -> ExtractionReport:
    with pytest.raises(ExtractionRefusedError) as info:
        extract_peptides(structure, chain, journal=Journal(), **kwargs)  # type: ignore[arg-type]
    assert isinstance(info.value.report, ExtractionReport)
    return info.value.report


# --- Jonctions ---------------------------------------------------------------------------


def test_continuous_chain_reports_each_junction_with_its_atoms_and_distance() -> None:
    journal = Journal()
    structure = Structure(atoms=_build((1, "ALA"), (2, "GLY"), (3, "SER")))
    result = extract_peptides(structure, "A", journal=journal)
    assert [r.code for r in result.peptides[0].residues] == ["ALA", "GLY", "SER"]
    first = result.report.junctions[0]
    assert (first.before, first.after) == ("A:ALA1", "A:GLY2")
    assert (first.measurements[0].c_atom, first.measurements[0].n_atom) == ("A:ALA1:C", "A:GLY2:N")
    assert first.distance == pytest.approx(BOND)
    assert first.category is JunctionCategory.CONTINUOUS
    assert first.decision is Decision.ACCEPT
    entry = journal.entries[-1]
    assert entry.action == "extract_peptide"
    assert entry.parameters["junctions"][0]["distance"] == pytest.approx(BOND)
    assert entry.parameters["criterion"]["max_c_n_distance"] == 1.8


@pytest.mark.parametrize(
    ("specs", "step"),
    [
        (((216, "ALA"), (218, "GLY")), 2),  # saut, comme ARG217 → GLY219 de 4ZHL
        (((1106, "ALA"), (219, "GLY")), -887),  # retour en arrière, comme 6WQA
        (((3, "GLY"), (3, "GLY", {"icode": "A"})), 0),  # code d'insertion
    ],
)
def test_numbering_discontinuity_with_a_bond_is_accepted_and_reported(
    specs: tuple[tuple[int, str] | tuple[int, str, dict[str, object]], ...], step: int
) -> None:
    result = extract_peptides(Structure(atoms=_build(*specs)), "A", journal=Journal())
    (junction,) = result.report.junctions
    assert junction.category is JunctionCategory.NUMBERING_DISCONTINUITY
    assert junction.decision is Decision.ACCEPT
    assert junction.numbering_step == step
    assert junction.distance == pytest.approx(BOND)
    assert len(result.peptides[0]) == 2


def test_undeclared_break_is_refused_with_its_distance() -> None:
    journal = Journal()
    structure = Structure(atoms=_build((1, "ALA"), (2, "GLY", {"gap": 14.1})))
    with pytest.raises(ExtractionRefusedError, match="break_undeclared") as info:
        extract_peptides(structure, "A", journal=journal)
    report = info.value.report
    assert isinstance(report, ExtractionReport)
    (junction,) = report.junctions
    assert junction.category is JunctionCategory.BREAK_UNDECLARED
    assert junction.decision is Decision.REJECT
    assert junction.distance == pytest.approx(14.1)
    assert journal.entries[-1].action == "extract_peptide_refused"
    assert journal.entries[-1].parameters["refusals"]


def test_declared_break_lists_the_missing_residues_between() -> None:
    structure = Structure(
        atoms=_build((1, "ALA"), (4, "GLY", {"gap": 9.0})),
        declared_missing_residues=(_rid(2, "SER"), _rid(3, "THR")),
    )
    report = _refused(structure)
    (junction,) = report.junctions
    assert junction.category is JunctionCategory.BREAK_DECLARED
    assert junction.declared_missing_between == ("A:SER2", "A:THR3")
    assert report.unobserved_termini == ()


def test_a_regular_numbering_does_not_hide_a_break() -> None:
    report = _refused(Structure(atoms=_build((1, "ALA"), (2, "GLY", {"gap": 5.0}))))
    assert report.junctions[0].category is JunctionCategory.BREAK_UNDECLARED


def test_missing_backbone_atom_makes_continuity_unmeasurable() -> None:
    structure = Structure(atoms=_build((1, "ALA", {"atoms": ("N", "CA", "O")}), (2, "GLY")))
    report = _refused(structure)
    junction = report.junctions[0]
    assert junction.category is JunctionCategory.CONTINUITY_UNMEASURABLE
    assert junction.measurements == ()
    assert "A:ALA1:C" in junction.evidence
    assert report.absent_backbone_atoms == ("A:ALA1:C",)


def test_discordant_alternate_conformations_make_continuity_unmeasurable() -> None:
    atoms = list(_build((1, "ALA", {"atoms": ("N", "CA", "O")}), (2, "GLY")))
    atoms += [
        Atom(90, "C", "C", "ALA", 1, "A", 2.0, 0.0, 0.0, alt_loc="A"),
        Atom(91, "C", "C", "ALA", 1, "A", -3.0, 0.0, 0.0, alt_loc="B"),
    ]
    report = _refused(Structure(atoms=tuple(atoms)))
    junction = report.junctions[0]
    assert junction.category is JunctionCategory.CONTINUITY_UNMEASURABLE
    assert len(junction.measurements) == 2
    assert {m.c_atom for m in junction.measurements} == {"A:ALA1:C[A]", "A:ALA1:C[B]"}


def test_threshold_can_be_imposed_and_the_distance_is_still_reported() -> None:
    structure = Structure(atoms=_build((1, "ALA"), (2, "GLY", {"gap": 2.0})))
    _refused(structure)
    criterion = load_peptide_bond_criterion(2.5)
    result = extract_peptides(structure, "A", journal=Journal(), criterion=criterion)
    assert result.report.junctions[0].distance == pytest.approx(2.0)
    assert "imposée" in result.report.criterion["source"]


# --- Résidus non observés et atomes absents ----------------------------------------------


def test_unobserved_terminal_residues_are_refused_by_default() -> None:
    structure = Structure(
        atoms=_build((2, "ALA"), (3, "GLY")),
        declared_missing_residues=(_rid(1, "MET"), _rid(4, "LYS")),
    )
    report = _refused(structure)
    assert [(t.side, t.residues, t.decision) for t in report.unobserved_termini] == [
        ("N", ("A:MET1",), Decision.REJECT),
        ("C", ("A:LYS4",), Decision.REJECT),
    ]
    assert report.segments[0].n_edge == "unobserved"
    assert report.segments[0].c_edge == "unobserved"


def test_missing_residue_of_another_chain_does_not_block() -> None:
    structure = Structure(
        atoms=_build((1, "ALA"), (2, "GLY")),
        declared_missing_residues=(_rid(1, "MET", chain="B"),),
    )
    assert len(structure_to_peptide(structure, "A", journal=Journal())) == 2


def test_declared_missing_residue_inside_a_bonded_chain_is_a_contradiction() -> None:
    structure = Structure(
        atoms=_build((1, "ALA"), (3, "GLY")),
        declared_missing_residues=(_rid(2, "SER"),),
    )
    report = _refused(structure)
    assert report.junctions[0].category is JunctionCategory.NUMBERING_DISCONTINUITY
    assert any("sans place" in reason for reason in report.refusals)


def test_declared_missing_atoms_and_absent_backbone_atoms_are_reported_separately() -> None:
    structure = Structure(
        atoms=_build((1, "ALA", {"atoms": ("N", "CA", "C")}), (2, "LYS")),
        declared_missing_atoms=(MissingAtom(_rid(2, "LYS"), "NZ"),),
    )
    report = analyze_chain(structure, "A")
    assert report.declared_missing_atoms == ("A:LYS2:NZ",)
    assert report.absent_backbone_atoms == ("A:ALA1:O",)
    assert report.accepted  # rapportés, pas bloquants : l'étape 4 les traite


# --- Composants --------------------------------------------------------------------------


def _ligand(name: str, number: int, atoms: int, chain: str = "A", x: float = 50.0) -> list[Atom]:
    return [
        Atom(500 + i, f"C{i + 1}", "C", name, number, chain, x + i, 10.0, 0.0, hetero=True)
        for i in range(atoms)
    ]


def test_unbound_ligand_sharing_the_chain_is_excluded_and_listed() -> None:
    atoms = (*_build((1, "ALA"), (2, "GLY")), *_ligand("ZMA", 1202, 3))
    result = extract_peptides(Structure(atoms=atoms), "A", journal=Journal())
    assert [r.code for r in result.peptides[0].residues] == ["ALA", "GLY"]
    (ligand,) = result.report.components
    assert ligand.category is ComponentCategory.UNBOUND_LIGAND
    assert ligand.decision is Decision.EXCLUDE
    assert ligand.min_distance_to_chain is not None and ligand.min_distance_to_chain > 10
    assert "A:ZMA1202" not in result.report.residues


def test_covalently_linked_glycans_are_side_chain_modifications_and_refused() -> None:
    atoms = (
        *_build((1, "ALA"), (2, "ASN", {"extra": ("CB", "CG", "ND2")})),
        *_ligand("NAG", 901, 14),
        *_ligand("NAG", 902, 14, x=70.0),
    )
    connections = (
        Connection(ConnectionKind.COVALENT, _rid(2, "ASN"), "ND2", _rid(901, "NAG"), "C1"),
        Connection(ConnectionKind.COVALENT, _rid(901, "NAG"), "O4", _rid(902, "NAG"), "C1"),
    )
    report = _refused(Structure(atoms=atoms, connections=connections))
    categories = {c.residue: c.category for c in report.components}
    assert categories == {
        "A:NAG901": ComponentCategory.SIDE_CHAIN_MODIFICATION,
        "A:NAG902": ComponentCategory.SIDE_CHAIN_MODIFICATION,  # lié par l'autre NAG
    }
    assert all(c.decision is Decision.REJECT for c in report.components)
    assert report.junctions[0].category is JunctionCategory.CONTINUOUS
    assert "A:NAG901" not in report.residues


def test_coordinated_ion_is_excluded_and_its_bonds_listed() -> None:
    atoms = (
        *_build((1, "HIS", {"extra": ("NE2",)}), (2, "GLY")),
        Atom(600, "ZN", "ZN", "ZN", 901, "A", 1.0, 5.0, 0.0, hetero=True),
    )
    metal = Connection(ConnectionKind.METAL, _rid(1, "HIS"), "NE2", _rid(901, "ZN"), "ZN")
    result = extract_peptides(Structure(atoms=atoms, connections=(metal,)), "A", journal=Journal())
    (ion,) = result.report.components
    assert ion.category is ComponentCategory.ION
    assert ion.decision is Decision.EXCLUDE
    assert ion.bonds == ("A:HIS1:NE2–A:ZN901:ZN (metal_coordination)",)
    assert len(result.peptides[0]) == 2


def test_modified_residue_in_the_chain_is_refused_with_its_position() -> None:
    report = _refused(Structure(atoms=_build((1, "ALA"), (2, "MSE"), (3, "GLY"))))
    (modified,) = report.components
    assert (modified.residue, modified.category) == ("A:MSE2", ComponentCategory.MODIFIED_RESIDUE)
    assert [j.category for j in report.junctions] == [JunctionCategory.CONTINUOUS] * 2


def test_terminal_cap_away_from_its_end_is_refused() -> None:
    atoms = _build((1, "ALA"), (2, "ACE", {"atoms": ("C", "O", "CH3")}), (3, "GLY"))
    report = _refused(Structure(atoms=atoms))
    assert any("groupe terminal N A:ACE2" in reason for reason in report.refusals)


def test_microheterogeneous_residue_is_refused() -> None:
    atoms = (
        *_build((1, "ALA")),
        Atom(9, "CA", "C", "SER", 2, "A", 0.0, 0.0, 0.0, alt_loc="A"),
        Atom(10, "CA", "C", "THR", 2, "A", 0.1, 0.0, 0.0, alt_loc="B"),
    )
    report = _refused(Structure(atoms=atoms))
    assert any("microhétérogénéité" in reason for reason in report.refusals)


# --- Liaisons déclarées ------------------------------------------------------------------


def test_interchain_disulfide_is_refused() -> None:
    atoms = (
        *_build((1, "CYS", {"extra": ("SG",)})),
        *_build((1, "CYS", {"extra": ("SG",)}), chain="B"),
    )
    link = Connection(
        ConnectionKind.DISULFIDE, _rid(1, "CYS"), "SG", _rid(1, "CYS", chain="B"), "SG"
    )
    report = _refused(Structure(atoms=atoms, connections=(link,)))
    assert any("autre chaîne" in reason for reason in report.refusals)


def test_declared_peptide_bond_is_accepted_but_other_covalent_links_are_refused() -> None:
    atoms = _build((1, "ALA"), (2, "GLY"), (3, "SER"))
    first, second, third = _rid(1, "ALA"), _rid(2, "GLY"), _rid(3, "SER")
    peptide_bond = Connection(ConnectionKind.COVALENT, first, "C", second, "N")
    structure = Structure(atoms=atoms, connections=(peptide_bond,))
    assert len(structure_to_peptide(structure, "A", journal=Journal())) == 3
    odd = Connection(ConnectionKind.COVALENT, first, "CA", third, "CA")
    report = _refused(Structure(atoms=atoms, connections=(odd,)))
    assert any("hors du squelette" in reason for reason in report.refusals)


# --- Fichiers d'essai écrits à la main ---------------------------------------------------


def _read(name: str, model: int | None = None) -> Structure:
    return read_structure(DATA / name, journal=Journal(), model=model)


@pytest.mark.parametrize("name", ["capped_peptide.pdb", "capped_peptide.cif"])
def test_capped_peptide_is_recognized_with_its_termini_and_disulfide(name: str) -> None:
    journal = Journal()
    result = extract_peptides(_read(name), "A", journal=journal)
    (peptide,) = result.peptides
    assert [r.code for r in peptide.residues] == ["ALA", "CYS", "GLY", "GLY", "CYS"]
    assert peptide.n_terminus.key == "acetyl"
    assert peptide.c_terminus.key == "amide"  # NH2 n'est pas pris pour un ion
    assert [(link.kind, link.positions) for link in peptide.links] == [(LinkKind.DISULFIDE, (2, 5))]
    assert N_TERMINUS_DEFAULTED not in peptide.metadata
    assert C_TERMINUS_DEFAULTED not in peptide.metadata
    assert peptide.metadata[STRUCTURE_CHAIN] == "A"
    assert peptide.metadata[STRUCTURE_SEGMENT] == "A:ACE0-A:NH25"
    categories = result.report.junctions_by_category()
    # Gly 3 → Gly 3A et Gly 3A → Cys 4 : code d'insertion de part et d'autre
    assert categories == {"continuous": 4, "numbering_discontinuity": 2}
    caps = [c.residue for c in result.report.components if c.category.value == "terminal_cap"]
    assert caps == ["A:ACE0", "A:NH25"]
    entry = journal.entries[-1]
    assert entry.parameters["peptides"][0]["n_terminus_from_structure"] is True
    assert entry.parameters["waters"] == 1


def test_free_termini_are_proven_by_h1_and_oxt() -> None:
    peptide = structure_to_peptide(_read("nmr_two_models.pdb", model=1), "B", journal=Journal())
    assert [r.code for r in peptide.residues] == ["GLY", "SER"]
    assert peptide.n_terminus.key == "free_amine"
    assert peptide.c_terminus.key == "free_acid"
    assert N_TERMINUS_DEFAULTED not in peptide.metadata
    assert C_TERMINUS_DEFAULTED not in peptide.metadata


def test_termini_without_evidence_are_defaulted_and_marked() -> None:
    peptide = structure_to_peptide(
        Structure(atoms=_build((1, "ALA"), (2, "GLY"))), "A", journal=Journal()
    )
    assert N_TERMINUS_DEFAULTED in peptide.metadata
    assert C_TERMINUS_DEFAULTED in peptide.metadata


def test_absent_chain_is_refused() -> None:
    with pytest.raises(StructureError, match="absente"):
        structure_to_peptide(_read("capped_peptide.pdb"), "Z", journal=Journal())


def test_waters_and_ions_are_not_part_of_the_peptide() -> None:
    result = extract_peptides(_read("capped_peptide.pdb"), "A", journal=Journal())
    assert len(result.peptides[0]) == 5
    excluded = {c.residue: c.category.value for c in result.report.components}
    assert excluded["A:HOH101"] == "water"
    assert excluded["A:ZN201"] == "ion"


# --- Décisions et rapport ----------------------------------------------------------------


def test_decisions_only_accept_their_own_values() -> None:
    assert ExtractionDecisions().to_dict() == {"breaks": "reject", "unobserved_termini": "reject"}
    with pytest.raises(ValueError, match="breaks"):
        ExtractionDecisions(breaks=Decision.OBSERVED_SEGMENT)
    with pytest.raises(ValueError, match="unobserved_termini"):
        ExtractionDecisions(unobserved_termini="split")  # type: ignore[arg-type]
    assert ExtractionDecisions(breaks="split").breaks is Decision.SPLIT  # type: ignore[arg-type]


def test_split_decision_is_recorded_on_each_break_and_cuts_segments() -> None:
    structure = Structure(atoms=_build((1, "ALA"), (2, "GLY", {"gap": 9.0}), (3, "SER")))
    report = analyze_chain(structure, "A", decisions=ExtractionDecisions(breaks=Decision.SPLIT))
    assert report.junctions[0].decision is Decision.SPLIT
    assert report.junctions[1].decision is Decision.ACCEPT
    assert [(s.first, s.last, s.n_edge, s.c_edge) for s in report.segments] == [
        ("A:ALA1", "A:ALA1", "chemical", "break"),
        ("A:GLY2", "A:SER3", "break", "chemical"),
    ]
    assert report.decisions["breaks"] == "split"
    assert not any("jonction" in reason for reason in report.refusals)


def test_report_is_serializable() -> None:
    structure = Structure(atoms=_build((1, "ALA"), (2, "GLY", {"gap": 9.0})))
    data = analyze_chain(structure, "A").to_dict()
    loaded = yaml.safe_load(yaml.safe_dump(data, allow_unicode=True))
    assert loaded["junctions"][0]["category"] == "break_undeclared"
    assert loaded["junctions"][0]["distance"] == 9.0
    assert loaded["accepted"] is False


# --- Structures réelles : 4ZHL, 6WQA, 7DDO -------------------------------------------------

REFERENCE_SHA256 = {
    "4ZHL.cif.gz": "b40655504ff4982cf8df49911019ada322d410a93d9bf66c65a9000e2531f60f",
    "6WQA.cif.gz": "484ff1689177aee680604631b9139199db8064f7de783e267f51cdfdf1769e6f",
    "7DDO.pdb.gz": "d6f4f7bacd3a8c8d21c4ec9d2543d74887a38fe469fd9840f8b9275f4039dd50",
}


@cache
def _reference(name: str) -> Structure:
    path = REFERENCE / name
    assert hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() == REFERENCE_SHA256[name]
    return read_structure(path, journal=Journal())


def _junction(report: ExtractionReport, before: str, after: str) -> object:
    (found,) = [j for j in report.junctions if (j.before, j.after) == (before, after)]
    return found


def test_4zhl_chymotrypsin_numbering_is_not_a_break() -> None:
    """Défaut 1 : le saut ARG217 → GLY219 était refusé comme « chaîne interrompue »."""
    result = extract_peptides(_reference("4ZHL.cif.gz"), "U", journal=Journal())
    report = result.report
    assert report.accepted
    assert len(result.peptides) == 1 and len(result.peptides[0]) == 247
    jump = _junction(report, "U:ARG217", "U:GLY219")
    assert jump.category is JunctionCategory.NUMBERING_DISCONTINUITY  # type: ignore[attr-defined]
    assert jump.numbering_step == 2  # type: ignore[attr-defined]
    assert jump.distance == pytest.approx(1.316, abs=1e-3)  # type: ignore[attr-defined]
    assert report.junctions_by_category() == {"continuous": 218, "numbering_discontinuity": 28}
    with_icode = [j for j in report.junctions if j.insertion_code]
    assert len(with_icode) == 27
    assert all(j.category is JunctionCategory.NUMBERING_DISCONTINUITY for j in with_icode)
    assert max(j.distance or 0 for j in report.junctions) < 1.8


def test_4zhl_peptide_chain_keeps_its_disulfide() -> None:
    result = extract_peptides(_reference("4ZHL.cif.gz"), "P", journal=Journal())
    (peptide,) = result.peptides
    assert len(peptide) == 10
    assert [link.positions for link in peptide.links] == [(1, 10)]
    assert result.report.junctions_by_category() == {"continuous": 9}


def test_6wqa_ligand_no_longer_blocks_and_the_undeclared_break_is_found() -> None:
    """Défaut 2 et lacune : `ZMA` comptait comme résidu, et la rupture n'était pas vue."""
    report = _refused(_reference("6WQA.cif.gz"))
    (zma,) = [c for c in report.components if c.name == "ZMA"]
    assert zma.category is ComponentCategory.UNBOUND_LIGAND
    assert zma.decision is Decision.EXCLUDE
    assert zma.min_distance_to_chain == pytest.approx(2.834, abs=1e-3)
    assert "A:ZMA1202" not in report.residues
    (reason,) = report.refusals  # seule la rupture bloque
    assert "A:ALA1043 → A:ASP1060" in reason
    gap = _junction(report, "A:ALA1043", "A:ASP1060")
    assert gap.category is JunctionCategory.BREAK_UNDECLARED  # type: ignore[attr-defined]
    assert gap.distance == pytest.approx(14.081, abs=1e-3)  # type: ignore[attr-defined]
    for before, after, step in (("A:LEU208", "A:ALA1001", 793), ("A:LEU1106", "A:GLU219", -887)):
        fusion = _junction(report, before, after)
        assert fusion.category is JunctionCategory.NUMBERING_DISCONTINUITY  # type: ignore[attr-defined]
        assert fusion.numbering_step == step  # type: ignore[attr-defined]
        assert fusion.distance < 1.8  # type: ignore[attr-defined]
    assert report.declared_missing_residues == ()


def test_6wqa_split_cuts_at_the_measured_break_only() -> None:
    report = analyze_chain(
        _reference("6WQA.cif.gz"), "A", decisions=ExtractionDecisions(breaks=Decision.SPLIT)
    )
    assert report.accepted
    assert [(s.first, s.last, s.residues) for s in report.segments] == [
        ("A:ASP-2", "A:ALA1043", 254),
        ("A:ASP1060", "A:LEU308", 137),
    ]


def test_7ddo_glycans_are_side_chain_modifications_not_residues() -> None:
    """Défaut 3 : les `NAG` liés aux Asn étaient comptés comme résidus de la chaîne."""
    report = _refused(_reference("7DDO.pdb.gz"))
    assert report.junctions_by_category() == {"continuous": 596}
    glycans = {c.residue: c.bonds for c in report.components if c.name == "NAG"}
    assert glycans == {
        "A:NAG902": ("A:ASN546:ND2–A:NAG902:C1 (covalent)",),
        "A:NAG903": ("A:ASN90:ND2–A:NAG903:C1 (covalent)",),
        "A:NAG904": ("A:ASN53:ND2–A:NAG904:C1 (covalent)",),
    }
    assert all(
        c.category is ComponentCategory.SIDE_CHAIN_MODIFICATION
        for c in report.components
        if c.name == "NAG"
    )
    assert not any(r.endswith(("NAG902", "NAG903", "NAG904")) for r in report.residues)
    assert len(report.refusals) == 3
    assert all("modification covalente" in reason for reason in report.refusals)
    (zinc,) = [c for c in report.components if c.name == "ZN"]
    assert zinc.category is ComponentCategory.ION
    assert zinc.decision is Decision.EXCLUDE
    assert len(zinc.bonds) == 3


def test_7ddo_chain_c_reports_unobserved_termini_where_they_are() -> None:
    report = _refused(_reference("7DDO.pdb.gz"), "C")
    sides = {t.side: t.residues for t in report.unobserved_termini}
    assert len(sides["N"]) == 14 and sides["N"][0] == "C:ARG319"
    assert sides["C"] == ("C:PRO527",)
    assert [c.residue for c in report.components] == ["C:NAG601"]
    assert report.junctions_by_category() == {"continuous": 193}


def test_alternate_conformations_are_paired_by_their_indicator() -> None:
    """Conformations A et B liées chacune ; les couples croisés A–B ne sont pas mesurés."""
    atoms = [
        *_build((1, "ALA", {"atoms": ("N", "CA", "O")}), (2, "GLY", {"atoms": ("CA", "C", "O")})),
        Atom(90, "C", "C", "ALA", 1, "A", 2.0, 0.0, 0.0, alt_loc="A"),
        Atom(91, "C", "C", "ALA", 1, "A", 2.0, 6.0, 0.0, alt_loc="B"),
        Atom(92, "N", "N", "GLY", 2, "A", 3.33, 0.0, 0.0, alt_loc="A"),
        Atom(93, "N", "N", "GLY", 2, "A", 3.33, 6.0, 0.0, alt_loc="B"),
    ]
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    (junction,) = report.junctions
    assert [(m.c_atom, m.n_atom) for m in junction.measurements] == [
        ("A:ALA1:C[A]", "A:GLY2:N[A]"),
        ("A:ALA1:C[B]", "A:GLY2:N[B]"),
    ]
    assert junction.category is JunctionCategory.CONTINUOUS


def test_covalent_link_to_another_chain_is_refused() -> None:
    atoms = (
        *_build((1, "LYS", {"extra": ("NZ",)}), (2, "GLY")),
        *_build((1, "GLY"), chain="B"),
    )
    link = Connection(ConnectionKind.COVALENT, _rid(1, "LYS"), "NZ", _rid(1, "GLY", chain="B"), "C")
    report = _refused(Structure(atoms=atoms, connections=(link,)))
    assert any("hors de la chaîne" in reason for reason in report.refusals)
