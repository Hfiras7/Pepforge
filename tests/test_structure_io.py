"""Tests de la lecture et de l'écriture des structures PDB et mmCIF (étape 3).

Fichiers d'essai écrits à la main dans `tests/data/structures/` :

- `capped_peptide.pdb` et `capped_peptide.cif` : la même molécule dans les deux
  formats, Ac-Ala-Cys-Gly-Gly(3A)-Cys-NH2, un pont disulfure 2–4 déclaré, une
  conformation alternative (CB d'Ala 1, A 0,60 / B 0,40), un code d'insertion
  (Gly 3A), une eau et un ion zinc ;
- `nmr_two_models.pdb` : ensemble de deux modèles de Gly-Ser, `NUMMDL 2`.

Les défauts sont injectés dans des copies, un par test.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from pepforge.core.structure import (
    STRUCTURE_FORMAT,
    STRUCTURE_MODEL_COUNT,
    STRUCTURE_SHA256,
    Atom,
    Connection,
    ConnectionKind,
    ResidueId,
    Structure,
)
from pepforge.exceptions import (
    RetrievalError,
    StructureError,
    StructureFormatError,
    UnsupportedFeatureError,
)
from pepforge.io import StructureFormat, read_structure, write_structure
from pepforge.io.structure_formats import GEMMI_VERSION, detect_format, format_from_name
from pepforge.report.journal import Journal

DATA = Path(__file__).parent / "data" / "structures"
PDB_FILE = DATA / "capped_peptide.pdb"
CIF_FILE = DATA / "capped_peptide.cif"
NMR_FILE = DATA / "nmr_two_models.pdb"


def _read(path: Path, **kwargs: object) -> Structure:
    return read_structure(path, journal=Journal(), **kwargs)  # type: ignore[arg-type]


def _copy(tmp_path: Path, source: Path, old: str, new: str, name: str | None = None) -> Path:
    text = source.read_text(encoding="utf-8")
    assert old in text, f"motif absent du fichier d'essai : {old!r}"
    target = tmp_path / (name or source.name)
    target.write_text(text.replace(old, new, 1), encoding="utf-8")
    return target


def _key(atom: Atom) -> tuple[object, ...]:
    return (
        atom.chain_id,
        atom.residue_number,
        atom.insertion_code,
        atom.residue_name,
        atom.name,
        atom.alt_loc,
        round(atom.x, 3),
        round(atom.y, 3),
        round(atom.z, 3),
        atom.element,
        atom.occupancy,
        atom.b_factor,
        atom.hetero,
    )


# --- Lecture correcte -------------------------------------------------------------------


def test_pdb_is_read_without_any_correction() -> None:
    structure = _read(PDB_FILE)
    assert len(structure) == 32
    assert structure.metadata[STRUCTURE_FORMAT] == "pdb"
    assert structure.metadata[STRUCTURE_MODEL_COUNT] == "1"
    assert structure.model is None
    # Conformation alternative gardée, avec ses deux occupations.
    cb = [a for a in structure.atoms if a.residue_number == 1 and a.name == "CB"]
    assert sorted((a.alt_loc, a.occupancy) for a in cb) == [("A", 0.6), ("B", 0.4)]
    # Code d'insertion gardé : Gly 3 et Gly 3A sont deux résidus distincts.
    glycines = {a.residue_id for a in structure.atoms if a.residue_name == "GLY"}
    assert glycines == {ResidueId("A", 3, "", "GLY"), ResidueId("A", 3, "A", "GLY")}
    # Eau, ion et groupes terminaux gardés, enregistrements HETATM signalés.
    names = {a.residue_name for a in structure.atoms if a.hetero}
    assert names == {"ACE", "NH2", "HOH", "ZN"}
    zinc = next(a for a in structure.atoms if a.residue_name == "ZN")
    assert zinc.element == "ZN"


def test_pdb_and_mmcif_of_the_same_molecule_give_the_same_atoms() -> None:
    pdb = _read(PDB_FILE)
    cif = _read(CIF_FILE)
    assert cif.metadata[STRUCTURE_FORMAT] == "mmcif"
    assert sorted(map(_key, pdb.atoms)) == sorted(map(_key, cif.atoms))
    assert pdb.connections == cif.connections


def test_declared_disulfide_is_kept() -> None:
    structure = _read(PDB_FILE)
    assert structure.connections == (
        Connection(
            ConnectionKind.DISULFIDE,
            ResidueId("A", 2, "", "CYS"),
            "SG",
            ResidueId("A", 4, "", "CYS"),
            "SG",
        ),
    )


def test_reading_is_journaled_with_its_provenance() -> None:
    journal = Journal()
    structure = read_structure(PDB_FILE, journal=journal)
    (entry,) = journal.entries
    assert entry.action == "read_structure"
    params = entry.parameters
    assert params["mode"] == "local"
    assert params["origin"] == str(PDB_FILE.resolve())
    assert params["size_bytes"] == PDB_FILE.stat().st_size
    assert len(params["sha256"]) == 64
    assert params["sha256"] == structure.metadata[STRUCTURE_SHA256]
    assert params["format"] == "pdb"
    assert params["reader"] == f"gemmi {GEMMI_VERSION}"
    assert params["atoms"] == 32


def test_gzip_file_is_read_and_its_received_bytes_are_hashed(tmp_path: Path) -> None:
    import hashlib

    raw = gzip.compress(PDB_FILE.read_bytes())
    target = tmp_path / "capped_peptide.pdb.gz"
    target.write_bytes(raw)
    journal = Journal()
    structure = read_structure(target, journal=journal)
    assert sorted(map(_key, structure.atoms)) == sorted(map(_key, _read(PDB_FILE).atoms))
    assert journal.entries[0].parameters["gzip"] is True
    assert structure.metadata[STRUCTURE_SHA256] == hashlib.sha256(raw).hexdigest()


# --- Modèles multiples ------------------------------------------------------------------


def test_multi_model_file_requires_an_explicit_model() -> None:
    with pytest.raises(StructureFormatError, match="2 modèles"):
        _read(NMR_FILE)


def test_chosen_model_is_read_and_recorded() -> None:
    first = _read(NMR_FILE, model=1)
    second = _read(NMR_FILE, model=2)
    assert first.model == 1 and second.model == 2
    assert first.metadata[STRUCTURE_MODEL_COUNT] == "2"
    shifts = {round(b.x - a.x, 3) for a, b in zip(first.atoms, second.atoms, strict=True)}
    assert shifts == {0.25}


def test_absent_model_is_refused() -> None:
    with pytest.raises(StructureFormatError, match="modèle 3 absent"):
        _read(NMR_FILE, model=3)


def test_truncated_ensemble_is_refused(tmp_path: Path) -> None:
    text = NMR_FILE.read_text(encoding="utf-8")
    truncated = text[: text.index("MODEL        2")] + "END\n"
    target = tmp_path / "truncated.pdb"
    target.write_text(truncated, encoding="utf-8")
    with pytest.raises(StructureFormatError, match="NUMMDL annonce 2"):
        _read(target)


def test_model_without_endmdl_is_refused(tmp_path: Path) -> None:
    target = _copy(tmp_path, NMR_FILE, "ENDMDL\nMODEL        2", "MODEL        2")
    with pytest.raises(StructureFormatError, match="MODEL sans ENDMDL"):
        _read(target, model=1)


# --- Lecture stricte : ce que gemmi laisserait passer ------------------------------------


def test_unreadable_pdb_coordinate_is_refused(tmp_path: Path) -> None:
    # gemmi 0.7.5 lit « 1x.104 » comme 1,0 sans erreur.
    target = _copy(tmp_path, PDB_FILE, "   0.800   2.500", "   0.800  2x.500")
    with pytest.raises(StructureFormatError, match="coordonnée y illisible"):
        _read(target)


def test_unreadable_mmcif_coordinate_is_refused(tmp_path: Path) -> None:
    # gemmi 0.7.5 lit « 1x.0 » et « ? » comme NaN.
    target = _copy(tmp_path, CIF_FILE, " 0.800 2.500 0.000 ", " 0.800 ? 0.000 ")
    with pytest.raises(StructureFormatError, match="coordonnée y illisible"):
        _read(target)


def test_missing_occupancy_is_left_unknown_not_filled(tmp_path: Path) -> None:
    # gemmi 0.7.5 remplit une occupation absente par 1,0 et un facteur B absent par 20,0.
    line = next(
        line for line in PDB_FILE.read_text().splitlines() if line.startswith("ATOM      5")
    )
    target = _copy(tmp_path, PDB_FILE, line, line[:54])
    journal = Journal()
    structure = read_structure(target, journal=journal)
    atom = next(a for a in structure.atoms if a.serial == 5)
    assert atom.occupancy is None and atom.b_factor is None
    assert journal.entries[0].parameters["atoms_without_occupancy"] == 1


def test_missing_element_is_counted(tmp_path: Path) -> None:
    line = next(
        line for line in PDB_FILE.read_text().splitlines() if line.startswith("ATOM      5")
    )
    target = _copy(tmp_path, PDB_FILE, line, line[:66])
    journal = Journal()
    read_structure(target, journal=journal)
    assert journal.entries[0].parameters["atoms_without_element"] == 1


def test_short_atom_line_is_refused(tmp_path: Path) -> None:
    line = next(
        line for line in PDB_FILE.read_text().splitlines() if line.startswith("ATOM      5")
    )
    target = _copy(tmp_path, PDB_FILE, line, line[:40])
    with pytest.raises(StructureFormatError, match="trop courte"):
        _read(target)


def test_duplicated_atom_is_refused(tmp_path: Path) -> None:
    line = next(
        line for line in PDB_FILE.read_text().splitlines() if line.startswith("ATOM      5")
    )
    target = _copy(tmp_path, PDB_FILE, line, line + "\n" + line)
    with pytest.raises(StructureFormatError, match="en double"):
        _read(target)


def test_file_without_atoms_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "empty.pdb"
    target.write_text("HEADER    NOTHING\nEND\n", encoding="utf-8")
    with pytest.raises(StructureFormatError, match="aucune ligne ATOM"):
        _read(target)


def test_mmcif_with_two_data_blocks_is_refused(tmp_path: Path) -> None:
    text = CIF_FILE.read_text(encoding="utf-8")
    target = tmp_path / "two.cif"
    target.write_text(text + text.replace("data_XTST", "data_OTHER"), encoding="utf-8")
    with pytest.raises(StructureFormatError, match="2 blocs"):
        _read(target)


def test_binary_content_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "binary.cif"
    target.write_bytes(b"\x89BCIF\xff\xfe\x00")
    with pytest.raises(StructureFormatError, match="non textuel"):
        _read(target)


# --- Format : extension, contenu, format demandé ----------------------------------------


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("x.pdb", StructureFormat.PDB),
        ("X.ENT", StructureFormat.PDB),
        ("x.cif", StructureFormat.MMCIF),
        ("x.mmcif.gz", StructureFormat.MMCIF),
        ("x.pdb.gz", StructureFormat.PDB),
        ("download", None),
        ("x.sdf", None),
    ],
)
def test_format_from_name(name: str, expected: StructureFormat | None) -> None:
    assert format_from_name(name) is expected


def test_format_is_recognized_from_content() -> None:
    assert detect_format(PDB_FILE.read_text()) is StructureFormat.PDB
    assert detect_format(CIF_FILE.read_text()) is StructureFormat.MMCIF
    assert detect_format("bonjour\n") is None


def test_extension_contradicting_content_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "wrong.cif"
    target.write_bytes(PDB_FILE.read_bytes())
    with pytest.raises(StructureFormatError, match="l'extension annonce mmcif"):
        _read(target)


def test_requested_format_contradicting_content_is_refused() -> None:
    with pytest.raises(StructureFormatError, match="le format demandé annonce mmcif"):
        _read(PDB_FILE, fmt="mmcif")


def test_unknown_requested_format_is_refused() -> None:
    with pytest.raises(StructureFormatError, match="Format inconnu"):
        _read(PDB_FILE, fmt="xyz")


def test_file_without_extension_is_read_from_its_content(tmp_path: Path) -> None:
    target = tmp_path / "download"
    target.write_bytes(CIF_FILE.read_bytes())
    journal = Journal()
    structure = read_structure(target, journal=journal)
    assert structure.metadata[STRUCTURE_FORMAT] == "mmcif"
    assert journal.entries[0].parameters["format_from_extension"] is None


def test_unrecognized_content_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "notes.txt"
    target.write_text("Ceci n'est pas une structure.\n", encoding="utf-8")
    with pytest.raises(StructureFormatError, match="contenu non reconnu"):
        _read(target)


def test_missing_local_file_is_refused(tmp_path: Path) -> None:
    with pytest.raises(RetrievalError, match="introuvable"):
        _read(tmp_path / "absent.pdb")


def test_structure_errors_share_a_base_class() -> None:
    assert issubclass(StructureFormatError, StructureError)
    assert issubclass(RetrievalError, StructureError)


# --- Résidus déclarés non observés ------------------------------------------------------

_REMARK_465 = """REMARK 465
REMARK 465 MISSING RESIDUES
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
REMARK 465
REMARK 465   M RES C SSSEQI
REMARK 465     MET A    -1
REMARK 465     SER A     6
"""


def test_remark_465_is_read(tmp_path: Path) -> None:
    target = _copy(tmp_path, PDB_FILE, "SSBOND", _REMARK_465 + "SSBOND")
    structure = _read(target)
    assert structure.declared_missing_residues == (
        ResidueId("A", -1, "", "MET"),
        ResidueId("A", 6, "", "SER"),
    )


def test_unreadable_remark_465_is_refused(tmp_path: Path) -> None:
    bad = _REMARK_465.replace("REMARK 465     SER A     6", "REMARK 465     SER A  six")
    target = _copy(tmp_path, PDB_FILE, "SSBOND", bad + "SSBOND")
    with pytest.raises(StructureFormatError, match="REMARK 465 illisible"):
        _read(target)


_UNOBS = """loop_
_pdbx_unobs_or_zero_occ_residues.id
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.polymer_flag
_pdbx_unobs_or_zero_occ_residues.occupancy_flag
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
_pdbx_unobs_or_zero_occ_residues.PDB_ins_code
1 1 Y 1 A MET -1 ?
2 1 Y 0 A ALA 1 ?
3 1 N 1 A HOH 300 ?
#
"""


def test_mmcif_unobserved_residues_are_read(tmp_path: Path) -> None:
    target = _copy(tmp_path, CIF_FILE, "loop_\n_atom_site", _UNOBS + "loop_\n_atom_site")
    structure = _read(target)
    # Seul le résidu de polymère non observé compte (occupancy_flag 0 : présent,
    # d'occupation nulle ; polymer_flag N : hors polymère).
    assert structure.declared_missing_residues == (ResidueId("A", -1, "", "MET"),)


_REMARK_470 = """REMARK 470
REMARK 470 MISSING ATOM
REMARK 470 THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;
REMARK 470 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;
REMARK 470 I=INSERTION CODE):
REMARK 470   M RES CSSEQI  ATOMS
REMARK 470     CYS A   2    SG   CB
REMARK 470     GLY A   3A   O
"""


def test_remark_470_missing_atoms_are_read(tmp_path: Path) -> None:
    target = _copy(tmp_path, PDB_FILE, "SSBOND", _REMARK_470 + "SSBOND")
    journal = Journal()
    structure = read_structure(target, journal=journal)
    assert [m.label() for m in structure.declared_missing_atoms] == [
        "A:CYS2:SG",
        "A:CYS2:CB",
        "A:GLY3A:O",
    ]
    assert journal.entries[0].parameters["declared_missing_atoms"] == 3


def test_unreadable_remark_470_is_refused(tmp_path: Path) -> None:
    bad = _REMARK_470.replace("REMARK 470     GLY A   3A   O", "REMARK 470     GLY A three O")
    target = _copy(tmp_path, PDB_FILE, "SSBOND", bad + "SSBOND")
    with pytest.raises(StructureFormatError, match="REMARK 470 illisible"):
        _read(target)


_UNOBS_ATOMS = """loop_
_pdbx_unobs_or_zero_occ_atoms.id
_pdbx_unobs_or_zero_occ_atoms.polymer_flag
_pdbx_unobs_or_zero_occ_atoms.occupancy_flag
_pdbx_unobs_or_zero_occ_atoms.PDB_model_num
_pdbx_unobs_or_zero_occ_atoms.auth_asym_id
_pdbx_unobs_or_zero_occ_atoms.auth_comp_id
_pdbx_unobs_or_zero_occ_atoms.auth_seq_id
_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code
_pdbx_unobs_or_zero_occ_atoms.auth_atom_id
1 Y 1 1 A CYS 2 ? SG
2 Y 0 1 A CYS 2 ? CB
#
"""


def test_mmcif_unobserved_atoms_are_read(tmp_path: Path) -> None:
    target = _copy(tmp_path, CIF_FILE, "loop_\n_atom_site", _UNOBS_ATOMS + "loop_\n_atom_site")
    structure = _read(target)
    # occupancy_flag 0 : atome présent d'occupation nulle, pas un atome manquant.
    assert [m.label() for m in structure.declared_missing_atoms] == ["A:CYS2:SG"]


# --- Écriture ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", ["pdb", "mmcif"])
@pytest.mark.parametrize("source", [PDB_FILE, CIF_FILE])
def test_write_then_read_back_loses_nothing(tmp_path: Path, source: Path, fmt: str) -> None:
    original = _read(source)
    journal = Journal()
    target = write_structure(original, tmp_path / f"out.{fmt}", fmt, journal=journal)
    again = _read(target)
    assert sorted(map(_key, again.atoms)) == sorted(map(_key, original.atoms))
    assert again.connections == original.connections
    (entry,) = journal.entries
    assert entry.action == "write_structure"
    assert entry.parameters["atoms"] == 32


def test_written_file_is_read_back_from_disk(tmp_path: Path) -> None:
    target = write_structure(_read(PDB_FILE), tmp_path / "out.pdb", "pdb", journal=Journal())
    assert target.read_text(encoding="utf-8").count("ATOM  ") > 0
    assert not list(tmp_path.glob(".*.tmp"))


def _atom(serial: int, chain: str = "A", number: int = 1) -> Atom:
    return Atom(serial, "CA", "C", "ALA", number, chain, 0.0, 0.0, float(serial))


def test_pdb_refuses_a_two_character_chain(tmp_path: Path) -> None:
    structure = Structure(atoms=(_atom(1, chain="AB"),))
    with pytest.raises(UnsupportedFeatureError, match="plus d'un caractère"):
        write_structure(structure, tmp_path / "out.pdb", "pdb", journal=Journal())
    assert not (tmp_path / "out.pdb").exists()


def test_mmcif_accepts_what_pdb_cannot_hold(tmp_path: Path) -> None:
    structure = Structure(atoms=(_atom(1, chain="AB", number=12345),))
    target = write_structure(structure, tmp_path / "out.cif", "mmcif", journal=Journal())
    (atom,) = _read(target).atoms
    assert (atom.chain_id, atom.residue_number) == ("AB", 12345)


def test_pdb_refuses_a_residue_number_out_of_range(tmp_path: Path) -> None:
    structure = Structure(atoms=(_atom(1, number=10000),))
    with pytest.raises(UnsupportedFeatureError, match="numéros de résidu"):
        write_structure(structure, tmp_path / "out.pdb", "pdb", journal=Journal())


def test_pdb_refuses_hydrogen_bonds_it_cannot_write(tmp_path: Path) -> None:
    rid = ResidueId("A", 1, "", "ALA")
    structure = Structure(
        atoms=(_atom(1),),
        connections=(Connection(ConnectionKind.HYDROGEN, rid, "CA", rid, "CA"),),
    )
    with pytest.raises(UnsupportedFeatureError, match="sans équivalent PDB"):
        write_structure(structure, tmp_path / "out.pdb", "pdb", journal=Journal())


def test_unknown_occupancy_is_written_as_the_format_default_and_journaled(
    tmp_path: Path,
) -> None:
    journal = Journal()
    structure = Structure(atoms=(_atom(1), _atom(2, number=2)))
    write_structure(structure, tmp_path / "out.pdb", "pdb", journal=journal)
    assert journal.entries[0].parameters["occupancy_written_as_default"] == 2
    assert {a.occupancy for a in _read(tmp_path / "out.pdb").atoms} == {1.0}


def test_empty_structure_is_not_written(tmp_path: Path) -> None:
    with pytest.raises(StructureFormatError, match="vide"):
        write_structure(Structure(atoms=()), tmp_path / "out.pdb", "pdb", journal=Journal())


def test_atom_refuses_non_finite_coordinates() -> None:
    with pytest.raises(ValueError, match="non finie"):
        Atom(1, "CA", "C", "ALA", 1, "A", float("nan"), 0.0, 0.0)


# --- Inventaire attendu contre lecture de gemmi : la comparaison elle-même ---------------


def _gemmi_and_raw(path: Path) -> tuple[object, list[object]]:
    import gemmi

    from pepforge.io.structure_formats import _scan_pdb

    text = path.read_text(encoding="utf-8")
    raw, _ = _scan_pdb(text, str(path))
    return gemmi.read_pdb_string(text), list(raw)


def test_atom_lost_by_the_reader_is_detected() -> None:
    import dataclasses

    from pepforge.io.structure_formats import _build_models

    structure, raw = _gemmi_and_raw(PDB_FILE)
    extra = dataclasses.replace(raw[0], key=(1, "A", 99, "", "ALA", "CA", ""))  # type: ignore[type-var]
    with pytest.raises(StructureFormatError, match="absent\\(s\\) de la lecture par gemmi"):
        _build_models(structure, [*raw, extra], "essai")  # type: ignore[arg-type, list-item]


def test_atom_invented_by_the_reader_is_detected() -> None:
    from pepforge.io.structure_formats import _build_models

    structure, raw = _gemmi_and_raw(PDB_FILE)
    with pytest.raises(StructureFormatError, match="atome absent du fichier"):
        _build_models(structure, raw[1:], "essai")  # type: ignore[arg-type]


def test_coordinate_changed_by_the_reader_is_detected() -> None:
    import dataclasses

    from pepforge.io.structure_formats import _build_models

    structure, raw = _gemmi_and_raw(PDB_FILE)
    shifted = dataclasses.replace(raw[0], x=raw[0].x + 0.01)  # type: ignore[attr-defined, type-var]
    with pytest.raises(StructureFormatError, match="coordonnées de gemmi"):
        _build_models(structure, [shifted, *raw[1:]], "essai")  # type: ignore[arg-type, list-item]


def test_writer_that_drops_an_atom_is_caught_and_nothing_is_written(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pepforge.io import structure_formats

    original = structure_formats.WRITERS[StructureFormat.PDB]

    def lossy(atoms: tuple[Atom, ...], connections: tuple[Connection, ...]) -> str:
        return original(atoms[:-1], connections)

    monkeypatch.setitem(structure_formats.WRITERS, StructureFormat.PDB, lossy)
    with pytest.raises(StructureFormatError, match="relus pour"):
        write_structure(_read(PDB_FILE), tmp_path / "out.pdb", "pdb", journal=Journal())
    assert not (tmp_path / "out.pdb").exists()
    assert not list(tmp_path.iterdir())


def test_atom_produced_twice_by_the_reader_is_detected() -> None:
    from pepforge.io.structure_formats import _build_models

    structure, raw = _gemmi_and_raw(PDB_FILE)
    residue = structure[0][0][0]  # type: ignore[index]
    residue.add_atom(residue[0])
    with pytest.raises(StructureFormatError, match="atome absent du fichier"):
        _build_models(structure, raw, "essai")  # type: ignore[arg-type]


def test_last_model_without_endmdl_is_refused(tmp_path: Path) -> None:
    text = NMR_FILE.read_text(encoding="utf-8")
    head, _, tail = text.rpartition("ENDMDL\n")
    target = tmp_path / "cut.pdb"
    target.write_text(head + tail, encoding="utf-8")
    with pytest.raises(StructureFormatError, match="dernier MODEL sans ENDMDL"):
        _read(target, model=1)
