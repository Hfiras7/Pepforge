"""Tests de l'extraction d'un peptide depuis une structure, et de son rapport.

Deux familles : des chaînes construites à la main, où chaque catégorie est
provoquée et sa géométrie connue, et trois structures réelles (4ZHL, 6WQA,
7DDO, `tests/data/structures/reference/`) qui ont révélé les défauts corrigés.
"""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import replace
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
from pepforge.exceptions import (
    ConfigurationError,
    ExtractionRefusedError,
    StructureError,
    UnsupportedFeatureError,
)
from pepforge.io import (
    BondEvidence,
    Certainty,
    ChainDecisions,
    Decision,
    ExtractionDecisions,
    ExtractionReport,
    analyze_chain,
    extract_chains,
    extract_peptides,
    read_structure,
    structure_to_peptide,
)
from pepforge.io.peptide_extraction import (
    STRUCTURE_CHAIN,
    STRUCTURE_EDGES,
    STRUCTURE_RESIDUES,
    STRUCTURE_SEGMENT,
    ComponentCategory,
    JunctionCategory,
)
from pepforge.report.journal import Journal
from pepforge.residues import assign_protonation, theoretical_net_charge
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
        # O à 1,23 Å de C, angle O=C–N de 123° : O(i)–N(i+1) = 2,25 Å, comme dans une
        # vraie liaison peptidique (avant l'étape C, O était à 1,79 Å du N suivant, ce
        # que la détection des contacts covalents signale à juste titre).
        positions = {
            "N": (x, 0.0),
            "CA": (x + 1.0, 0.5),
            "C": (x + 2.0, 0.0),
            "O": (x + 1.33, 1.03),
        }
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


# Contacts SG–SG non déclarés trouvés à l'étape C dans ces copies de fichiers (le module 3
# extrayait ces chaînes sans les voir). Revue de Firas du 28/11/2026 : ils restent des
# contacts **possibles**, jamais admis comme ponts disulfure dans les tests (l'entrée
# originale de 6WQA les déclare, ce que ces copies sans en-tête ne peuvent pas dire) ; la chaîne U
# de 4ZHL reste refusée (la validation du fichier les classe en « close contact »).
ZHL_U_CONTACTS = ("U:CYS50:SG–U:CYS111:SG", "U:CYS136:SG–U:CYS201:SG")
WQA_A_CONTACTS = (
    "A:CYS71:SG–A:CYS159:SG",
    "A:CYS74:SG–A:CYS146:SG",
    "A:CYS77:SG–A:CYS166:SG",
    "A:CYS259:SG–A:CYS262:SG",
)


def _contact_refusals(report: ExtractionReport, contacts: tuple[str, ...]) -> list[str]:
    """Les refus de contacts possibles du rapport, un par contact attendu, dans l'ordre."""
    return [
        next(r for r in report.refusals if contact in r and "contact covalent possible" in r)
        for contact in contacts
    ]


def _without_sg(structure: Structure, residues: tuple[str, ...]) -> Structure:
    """Copie d'essai sans l'atome SG des résidus nommés.

    Sert seulement à isoler un autre mécanisme du module 3 (numérotation, rupture,
    décisions par chaîne) des contacts SG–SG non tranchés : aucun pont n'est admis,
    l'atome retiré est rapporté comme absent. Ce n'est jamais une préparation.
    """
    removed = [a for a in structure.atoms if a.name == "SG" and a.residue_id.label() in residues]
    assert len(removed) == len(residues)
    kept = tuple(a for a in structure.atoms if a not in removed)
    return replace(structure, atoms=kept)


ZHL_U_SG_REMOVED = ("U:CYS111", "U:CYS201")
WQA_A_SG_REMOVED = ("A:CYS159", "A:CYS146", "A:CYS166", "A:CYS262")


def _junction(report: ExtractionReport, before: str, after: str) -> object:
    (found,) = [j for j in report.junctions if (j.before, j.after) == (before, after)]
    return found


def test_4zhl_chymotrypsin_numbering_is_not_a_break() -> None:
    """Défaut 1 : le saut ARG217 → GLY219 était refusé comme « chaîne interrompue »."""
    report = _refused(_reference("4ZHL.cif.gz"), "U")
    # Seuls les deux contacts SG–SG non tranchés bloquent (revue de Firas, 28/11/2026).
    assert sorted(report.refusals) == sorted(_contact_refusals(report, ZHL_U_CONTACTS))
    assert len(report.residues) == 247
    jump = _junction(report, "U:ARG217", "U:GLY219")
    assert jump.category is JunctionCategory.NUMBERING_DISCONTINUITY  # type: ignore[attr-defined]
    assert jump.numbering_step == 2  # type: ignore[attr-defined]
    assert jump.distance == pytest.approx(1.316, abs=1e-3)  # type: ignore[attr-defined]
    assert report.junctions_by_category() == {"continuous": 218, "numbering_discontinuity": 28}
    with_icode = [j for j in report.junctions if j.insertion_code]
    assert len(with_icode) == 27
    assert all(j.category is JunctionCategory.NUMBERING_DISCONTINUITY for j in with_icode)
    assert max(j.distance or 0 for j in report.junctions) < 1.8


def test_4zhl_chain_u_stays_refused_and_its_contacts_stay_possible() -> None:
    """Revue de Firas : refus maintenu, les deux contacts ne sont pas des ponts."""
    report = _refused(_reference("4ZHL.cif.gz"), "U")
    contacts = {b.label(): b for b in report.bonds if b.label() in ZHL_U_CONTACTS}
    assert set(contacts) == set(ZHL_U_CONTACTS)
    assert {label: round(b.distance or 0, 3) for label, b in contacts.items()} == {
        "U:CYS50:SG–U:CYS111:SG": 1.637,
        "U:CYS136:SG–U:CYS201:SG": 1.515,
    }
    for bond in contacts.values():
        assert bond.evidence is BondEvidence.GEOMETRIC
        assert bond.certainty is Certainty.POSSIBLE


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
    # La rupture, puis les quatre contacts SG–SG laissés possibles (revue du 28/11/2026).
    contacts = _contact_refusals(report, WQA_A_CONTACTS)
    (reason,) = [r for r in report.refusals if r not in contacts]
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
    assert sorted(report.refusals) == sorted(_contact_refusals(report, WQA_A_CONTACTS))
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


# --- Décisions split et observed_segment : extrémité non chimique ----------------------


def test_observed_segment_gives_non_chemical_termini_where_residues_are_missing() -> None:
    journal = Journal()
    structure = Structure(
        atoms=_build((2, "ALA"), (3, "GLY", {"extra": ("OXT",)})),
        declared_missing_residues=(_rid(1, "MET"), _rid(4, "LYS")),
    )
    decisions = ExtractionDecisions(unobserved_termini=Decision.OBSERVED_SEGMENT)
    result = extract_peptides(structure, "A", journal=journal, decisions=decisions)
    (peptide,) = result.peptides
    assert peptide.n_terminus.key == peptide.c_terminus.key == "unobserved"
    assert not peptide.n_terminus.chemical and not peptide.c_terminus.chemical
    # Un OXT au bord d'un segment ne fait pas un acide libre : la chaîne continue.
    assert peptide.metadata[STRUCTURE_EDGES] == "unobserved/unobserved"
    assert N_TERMINUS_DEFAULTED not in peptide.metadata
    assert C_TERMINUS_DEFAULTED not in peptide.metadata
    assert [t.decision for t in result.report.unobserved_termini] == [Decision.OBSERVED_SEGMENT] * 2
    entry = journal.entries[-1]
    assert entry.parameters["decisions"]["unobserved_termini"] == "observed_segment"
    assert entry.parameters["peptides"][0]["edges"] == "unobserved/unobserved"
    assert entry.parameters["peptides"][0]["c_terminus_from_structure"] is False


def test_observed_segment_keeps_a_chemical_end_when_only_one_side_is_missing() -> None:
    structure = Structure(
        atoms=_build((2, "ALA"), (3, "GLY", {"extra": ("OXT",)})),
        declared_missing_residues=(_rid(1, "MET"),),
    )
    decisions = ExtractionDecisions(unobserved_termini=Decision.OBSERVED_SEGMENT)
    (peptide,) = extract_peptides(structure, "A", journal=Journal(), decisions=decisions).peptides
    assert peptide.n_terminus.key == "unobserved"
    assert peptide.c_terminus.key == "free_acid"  # OXT, vraie extrémité


def test_split_gives_one_peptide_per_segment_with_non_chemical_inner_edges() -> None:
    structure = Structure(
        atoms=_build((1, "ALA"), (2, "GLY"), (5, "SER", {"gap": 12.0}), (6, "THR"))
    )
    result = extract_peptides(
        structure, "A", journal=Journal(), decisions=ExtractionDecisions(breaks=Decision.SPLIT)
    )
    first, second = result.peptides
    assert [r.code for r in first.residues] == ["ALA", "GLY"]
    assert [r.code for r in second.residues] == ["SER", "THR"]
    assert first.c_terminus.key == "unobserved" and second.n_terminus.key == "unobserved"
    assert N_TERMINUS_DEFAULTED in first.metadata  # vraie extrémité, sans preuve : marquée
    assert C_TERMINUS_DEFAULTED in second.metadata
    assert first.metadata[STRUCTURE_SEGMENT] == "A:ALA1-A:GLY2"
    assert second.metadata[STRUCTURE_EDGES] == "break/chemical"


def test_split_refuses_a_disulfide_across_two_segments() -> None:
    atoms = _build(
        (1, "CYS", {"extra": ("SG",)}), (2, "GLY"), (3, "CYS", {"gap": 12.0, "extra": ("SG",)})
    )
    bridge = Connection(ConnectionKind.DISULFIDE, _rid(1, "CYS"), "SG", _rid(3, "CYS"), "SG")
    report = _refused(
        Structure(atoms=atoms, connections=(bridge,)),
        decisions=ExtractionDecisions(breaks=Decision.SPLIT),
    )
    assert any("entre deux segments" in reason for reason in report.refusals)


def test_6wqa_split_gives_two_peptides_at_the_undeclared_break() -> None:
    result = extract_peptides(
        _without_sg(_reference("6WQA.cif.gz"), WQA_A_SG_REMOVED),
        "A",
        journal=Journal(),
        decisions=ExtractionDecisions(breaks=Decision.SPLIT),
    )
    first, second = result.peptides
    assert (len(first), len(second)) == (254, 137)
    assert first.c_terminus.key == "unobserved" and second.n_terminus.key == "unobserved"
    assert "ZMA" not in {r.code for p in result.peptides for r in p.residues}


def test_7ddo_observed_segment_does_not_hide_the_glycan() -> None:
    report = _refused(
        _reference("7DDO.pdb.gz"),
        "C",
        decisions=ExtractionDecisions(unobserved_termini=Decision.OBSERVED_SEGMENT),
    )
    assert len(report.refusals) == 1
    assert "C:NAG601" in report.refusals[0]


# --- Validation finale : numérotation d'origine, décisions par chaîne et par catégorie ----


def test_4zhl_peptide_keeps_the_original_numbering_of_each_position() -> None:
    (peptide,) = extract_peptides(
        _without_sg(_reference("4ZHL.cif.gz"), ZHL_U_SG_REMOVED), "U", journal=Journal()
    ).peptides
    original = peptide.metadata[STRUCTURE_RESIDUES].split(",")
    assert len(original) == len(peptide) == 247
    position = original.index("U:ARG217")
    assert original[position + 1] == "U:GLY219"  # le saut reste visible, rien n'est renuméroté
    assert peptide.residues[position].code == "ARG"
    assert peptide.residues[position + 1].code == "GLY"
    assert any(label.endswith("A") and label[-2].isdigit() for label in original)  # 184A...


def test_numbering_discontinuity_evidence_names_the_jump_and_the_distance() -> None:
    report = analyze_chain(_reference("4ZHL.cif.gz"), "U")
    jump = _junction(report, "U:ARG217", "U:GLY219")
    assert "numérotation +2" in jump.evidence  # type: ignore[attr-defined]
    assert "1.316 Å ≤ 1.8 Å" in jump.evidence  # type: ignore[attr-defined]
    (measure,) = jump.measurements  # type: ignore[attr-defined]
    assert (measure.c_atom, measure.n_atom) == ("U:ARG217:C", "U:GLY219:N")


def test_conformations_without_a_compatible_pair_are_unmeasurable_not_bonded() -> None:
    """Une absence de mesure n'est jamais une preuve de continuité."""
    atoms = [
        *_build((1, "ALA", {"atoms": ("N", "CA", "O")}), (2, "GLY", {"atoms": ("CA", "C", "O")})),
        Atom(90, "C", "C", "ALA", 1, "A", 2.0, 0.0, 0.0, alt_loc="A"),
        Atom(91, "N", "N", "GLY", 2, "A", 3.33, 0.0, 0.0, alt_loc="B"),
    ]
    report = _refused(Structure(atoms=tuple(atoms)))
    (junction,) = report.junctions
    assert junction.category is JunctionCategory.CONTINUITY_UNMEASURABLE
    assert junction.measurements == ()
    assert "aucune paire mesurable" in junction.evidence
    assert junction.decision is Decision.REJECT


def _two_broken_chains() -> Structure:
    return Structure(
        atoms=(
            *_build((1, "ALA"), (2, "GLY"), (3, "SER", {"gap": 10.0})),
            # Chaîne B à 20 Å de la chaîne A : avant l'étape C, les deux chaînes étaient
            # superposées, ce que la détection des contacts entre chaînes signale.
            *(
                replace(atom, z=atom.z + 20.0)
                for atom in _build((1, "LYS"), (2, "GLY"), (3, "THR", {"gap": 10.0}), chain="B")
            ),
        )
    )


def test_decisions_of_one_chain_do_not_propagate_to_another() -> None:
    journal = Journal()
    decisions = ChainDecisions({"A": ExtractionDecisions(breaks=Decision.SPLIT)})
    assert decisions.for_chain("B") == ExtractionDecisions()
    with pytest.raises(ExtractionRefusedError, match="chaîne B") as info:
        extract_chains(_two_broken_chains(), ["A", "B"], journal=journal, decisions=decisions)
    report_b = info.value.report
    assert isinstance(report_b, ExtractionReport)
    assert report_b.decisions == {"breaks": "reject", "unobserved_termini": "reject"}
    assert report_b.junctions[1].decision is Decision.REJECT
    accepted, refused = journal.entries[-2], journal.entries[-1]
    assert (accepted.action, accepted.parameters["chain"]) == ("extract_peptide", "A")
    assert accepted.parameters["decisions"]["breaks"] == "split"
    assert len(accepted.parameters["peptides"]) == 2
    assert (refused.action, refused.parameters["chain"]) == ("extract_peptide_refused", "B")


def test_each_chain_receives_its_own_decision() -> None:
    decisions = ChainDecisions.from_dict({"A": {"breaks": "split"}, "B": {"breaks": "split"}})
    results = extract_chains(
        _two_broken_chains(), ["A", "B"], journal=Journal(), decisions=decisions
    )
    assert {chain: len(r.peptides) for chain, r in results.items()} == {"A": 2, "B": 2}
    assert results["B"].peptides[1].metadata[STRUCTURE_CHAIN] == "B"
    assert decisions.to_dict() == {
        "A": {"breaks": "split", "unobserved_termini": "reject"},
        "B": {"breaks": "split", "unobserved_termini": "reject"},
    }


def test_chain_decisions_refuse_unknown_chains_categories_and_values() -> None:
    structure = _two_broken_chains()
    stray = ChainDecisions({"Z": ExtractionDecisions(breaks=Decision.SPLIT)})
    with pytest.raises(ConfigurationError, match="non extraites"):
        extract_chains(structure, ["A"], journal=Journal(), decisions=stray)
    with pytest.raises(StructureError, match="absente"):
        extract_chains(structure, ["A", "Q"], journal=Journal())
    with pytest.raises(ConfigurationError, match="inconnue"):
        ChainDecisions.from_dict({"A": {"missing_residues": "split"}})
    with pytest.raises(ConfigurationError, match="chaîne A"):
        ChainDecisions.from_dict({"A": {"breaks": "observed_segment"}})


def _break_and_unobserved_ends() -> Structure:
    return Structure(
        atoms=_build((2, "ALA"), (3, "GLY", {"extra": ("OXT",)}), (4, "SER", {"gap": 10.0})),
        declared_missing_residues=(_rid(1, "MET"),),
    )


def test_a_decision_on_breaks_does_not_accept_unobserved_termini() -> None:
    report = _refused(
        _break_and_unobserved_ends(), decisions=ExtractionDecisions(breaks=Decision.SPLIT)
    )
    assert report.junctions[1].decision is Decision.SPLIT
    assert [t.decision for t in report.unobserved_termini] == [Decision.REJECT]
    (reason,) = report.refusals
    assert "extrémité N" in reason


def test_a_decision_on_unobserved_termini_does_not_accept_breaks() -> None:
    report = _refused(
        _break_and_unobserved_ends(),
        decisions=ExtractionDecisions(unobserved_termini=Decision.OBSERVED_SEGMENT),
    )
    assert [t.decision for t in report.unobserved_termini] == [Decision.OBSERVED_SEGMENT]
    assert report.junctions[1].decision is Decision.REJECT
    (reason,) = report.refusals
    assert "break_undeclared" in reason


def test_both_decisions_together_and_oxt_before_a_break_is_not_a_free_acid() -> None:
    decisions = ExtractionDecisions(
        breaks=Decision.SPLIT, unobserved_termini=Decision.OBSERVED_SEGMENT
    )
    first, second = extract_peptides(
        _break_and_unobserved_ends(), "A", journal=Journal(), decisions=decisions
    ).peptides
    assert (first.n_terminus.key, first.c_terminus.key) == ("unobserved", "unobserved")
    assert first.metadata[STRUCTURE_EDGES] == "unobserved/break"
    assert first.metadata[STRUCTURE_RESIDUES] == "A:ALA2,A:GLY3"
    assert (second.n_terminus.key, second.c_terminus.key) == ("unobserved", "free_acid")
    assert C_TERMINUS_DEFAULTED in second.metadata  # vraie extrémité sans preuve : marquée


def test_extracted_segment_cannot_be_protonated_or_charged() -> None:
    structure = Structure(
        atoms=_build((2, "LYS"), (3, "GLY")),
        declared_missing_residues=(_rid(1, "MET"),),
    )
    decisions = ExtractionDecisions(unobserved_termini=Decision.OBSERVED_SEGMENT)
    (peptide,) = extract_peptides(structure, "A", journal=Journal(), decisions=decisions).peptides
    with pytest.raises(UnsupportedFeatureError, match="pas une extrémité chimique"):
        assign_protonation(peptide, 7.4)
    with pytest.raises(UnsupportedFeatureError, match="coiffer ou reconstruire"):
        theoretical_net_charge(peptide, 7.4)


def test_4zhl_two_chains_are_extracted_independently() -> None:
    results = extract_chains(
        _without_sg(_reference("4ZHL.cif.gz"), ZHL_U_SG_REMOVED),
        ["U", "P"],
        journal=Journal(),
        decisions=ChainDecisions({"U": ExtractionDecisions()}),
    )
    assert {c: len(r.peptides[0]) for c, r in results.items()} == {"U": 247, "P": 10}
    assert results["P"].report.junctions_by_category() == {"continuous": 9}
