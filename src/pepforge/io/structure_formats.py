"""Formats de fichiers de structure : registre, détection, lecture stricte.

Architecture générique (demande de Firas, 27 novembre 2026) : chaque format
s'inscrit dans un registre (`READERS`, `WRITERS`) avec son lecteur et son
écrivain. Ajouter un format, c'est ajouter une valeur à `StructureFormat`, ses
extensions, sa signature de contenu, puis inscrire ses fonctions : rien d'autre
ne change, ni dans `io.structure`, ni ailleurs.

Bibliothèque de lecture : gemmi (Wojdyr, *J. Open Source Softw.* 7(73), 4200,
2022, DOI 10.21105/joss.04200), choisie pour garder conformations alternatives,
codes d'insertion et modèles multiples. Elle n'est importée **que dans ce
fichier** (vérifié par `tests/test_architecture.py`) : le reste de PepForge ne
voit que `Structure`.

Lecture stricte : gemmi tolère des fichiers qu'il faudrait refuser. Constaté sur
gemmi 0.7.5 (27 novembre 2026) :

- en PDB, la coordonnée `1x.104` est lue `1.0` sans erreur ;
- en PDB, une occupation absente devient 1,0 et un facteur B absent 20,0 ;
- en mmCIF, `1x.0` et `?` deviennent NaN ;
- un élément absent (colonnes 77-78 en PDB, `type_symbol` en mmCIF) est déduit
  du nom de l'atome.

PepForge relit donc lui-même chaque enregistrement d'atome du fichier brut (les
lignes `ATOM` / `HETATM` du PDB, la table `_atom_site` du mmCIF), avec des
règles strictes, et en tire un **inventaire attendu**. La structure produite
par gemmi est ensuite comparée à cet inventaire, atome par atome : nombre,
identité (modèle, chaîne, numéro, code d'insertion, résidu, nom, conformation)
et coordonnées. Le moindre écart lève `StructureFormatError` (méthode de la
section 17 de l'architecture : inventaire attendu contre fichier produit).
"""

from __future__ import annotations

import importlib.metadata
import math
import re
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from pathlib import PurePath

import gemmi

from pepforge.core.structure import Atom, Connection, ConnectionKind, ResidueId
from pepforge.exceptions import StructureFormatError

GEMMI_VERSION: str = importlib.metadata.version("gemmi")


class StructureFormat(StrEnum):
    """Formats de structure pris en charge."""

    PDB = "pdb"
    MMCIF = "mmcif"


EXTENSIONS: dict[str, StructureFormat] = {
    ".pdb": StructureFormat.PDB,
    ".ent": StructureFormat.PDB,
    ".cif": StructureFormat.MMCIF,
    ".mmcif": StructureFormat.MMCIF,
}
"""Extensions reconnues (comparées sans tenir compte de la casse, `.gz` retiré)."""

PDB_RECORDS = frozenset(
    {
        "HEADER",
        "TITLE ",
        "COMPND",
        "SOURCE",
        "REMARK",
        "CRYST1",
        "MODEL ",
        "ATOM  ",
        "HETATM",
        "SEQRES",
        "EXPDTA",
        "AUTHOR",
        "KEYWDS",
        "NUMMDL",
        "HELIX ",
        "SHEET ",
        "SSBOND",
        "LINK  ",
        "HETNAM",
        "HET   ",
        "DBREF ",
        "ORIGX1",
        "SCALE1",
    }
)
"""Noms d'enregistrements dont la présence en début de ligne signe un fichier PDB."""


def format_from_name(name: str) -> StructureFormat | None:
    """Format annoncé par l'extension du nom de fichier ; `None` si elle n'en annonce pas."""
    path = PurePath(name.lower())
    if path.suffix == ".gz":
        path = path.with_suffix("")
    return EXTENSIONS.get(path.suffix)


def decode_text(content: bytes, origin: str) -> str:
    """Décode le contenu : UTF-8 (dont l'ASCII des fichiers PDB), sinon refus."""
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise StructureFormatError(
            f"{origin} : contenu non textuel ou non UTF-8 (octet {exc.start}) ; "
            "un fichier binaire (BinaryCIF, MMTF, archive) n'est pas pris en charge"
        ) from None
    return text.lstrip("﻿")


def detect_format(text: str) -> StructureFormat | None:
    """Format reconnu à la signature du contenu ; `None` si aucune ne correspond."""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.lower().startswith("data_"):
            return StructureFormat.MMCIF
        break
    for line in text.splitlines()[:2000]:
        if line[:6].ljust(6) in PDB_RECORDS:
            return StructureFormat.PDB
    return None


# --- Inventaire attendu, tiré du fichier brut --------------------------------------------

_FLOAT = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
_INT = re.compile(r"[+-]?\d+")

AtomKey = tuple[int, str, int, str, str, str, str]
"""Identité d'un atome : modèle, chaîne, numéro, code d'insertion, résidu, nom, conformation."""


@dataclass(frozen=True)
class RawAtom:
    """Atome tel qu'écrit dans le fichier, relu par PepForge."""

    key: AtomKey
    x: float
    y: float
    z: float
    occupancy: float | None
    b_factor: float | None
    element: str | None
    hetero: bool


@dataclass(frozen=True)
class ParsedFile:
    """Ce qu'un lecteur de format rend : atomes vérifiés, par modèle, et déclarations du fichier."""

    format: StructureFormat
    models: dict[int, tuple[Atom, ...]]
    connections: tuple[Connection, ...]
    declared_missing: dict[int | None, tuple[ResidueId, ...]]
    """Résidus déclarés non observés, par modèle ; clé `None` : valables pour tous les modèles."""

    declared_model_count: int | None
    atoms_without_element: int
    """Atomes dont le fichier ne donne pas l'élément (déduit par gemmi du nom d'atome)."""

    atoms_without_occupancy: int
    """Atomes dont le fichier ne donne pas l'occupation (laissée à `None`)."""


def _strict_float(text: str, what: str, where: str) -> float:
    value = text.strip()
    if not _FLOAT.fullmatch(value):
        raise StructureFormatError(f"{where} : {what} illisible : {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise StructureFormatError(f"{where} : {what} non finie : {value!r}")
    return number


def _optional_float(text: str, what: str, where: str) -> float | None:
    value = text.strip()
    if value in ("", "?", "."):
        return None
    return _strict_float(value, what, where)


def _strict_int(text: str, what: str, where: str) -> int:
    value = text.strip()
    if not _INT.fullmatch(value):
        raise StructureFormatError(f"{where} : {what} illisible : {value!r}")
    return int(value)


def _scan_pdb(text: str, origin: str) -> tuple[list[RawAtom], int | None]:
    """Relit les lignes ATOM / HETATM selon les colonnes de la spécification PDB 3.3."""
    atoms: list[RawAtom] = []
    model = 1
    in_model = False
    declared: int | None = None
    for number, line in enumerate(text.splitlines(), start=1):
        record = line[:6]
        where = f"{origin}, ligne {number}"
        if record == "NUMMDL":
            declared = _strict_int(line[10:14], "nombre de modèles (NUMMDL)", where)
        elif record == "MODEL ":
            if in_model:
                raise StructureFormatError(f"{where} : MODEL sans ENDMDL pour le modèle précédent")
            model = _strict_int(line[10:14] if len(line) > 10 else "", "numéro de modèle", where)
            in_model = True
        elif record == "ENDMDL":
            in_model = False
        elif record in ("ATOM  ", "HETATM"):
            if len(line.rstrip()) < 54:
                raise StructureFormatError(
                    f"{where} : ligne d'atome trop courte ({len(line.rstrip())} colonnes, "
                    "54 au moins pour les coordonnées)"
                )
            padded = line.ljust(80)
            name = padded[12:16].strip()
            if not name:
                raise StructureFormatError(f"{where} : nom d'atome vide")
            residue_number = _strict_int(padded[22:26], "numéro de résidu", where)
            x = _strict_float(padded[30:38], "coordonnée x", where)
            y = _strict_float(padded[38:46], "coordonnée y", where)
            z = _strict_float(padded[46:54], "coordonnée z", where)
            occupancy = _optional_float(padded[54:60], "occupation", where)
            b_factor = _optional_float(padded[60:66], "facteur B", where)
            element = padded[76:78].strip().upper() or None
            key: AtomKey = (
                model,
                padded[21].strip(),
                residue_number,
                padded[26].strip(),
                padded[17:20].strip(),
                name,
                padded[16].strip(),
            )
            atoms.append(RawAtom(key, x, y, z, occupancy, b_factor, element, record == "HETATM"))
    if in_model:
        raise StructureFormatError(f"{origin} : dernier MODEL sans ENDMDL (fichier tronqué ?)")
    return atoms, declared


_ATOM_SITE_REQUIRED = ("Cartn_x", "Cartn_y", "Cartn_z")


def _scan_mmcif(block: gemmi.cif.Block, origin: str) -> list[RawAtom]:
    """Relit la table `_atom_site` du bloc mmCIF, valeur par valeur."""
    table = block.find_mmcif_category("_atom_site.")
    if len(table) == 0:
        raise StructureFormatError(f"{origin} : aucune table _atom_site (aucun atome)")
    tags = {tag.removeprefix("_atom_site."): i for i, tag in enumerate(table.tags)}
    missing = [t for t in _ATOM_SITE_REQUIRED if t not in tags]
    if missing:
        raise StructureFormatError(f"{origin} : colonnes _atom_site absentes : {missing}")

    def col(*names: str) -> int | None:
        for name in names:
            if name in tags:
                return tags[name]
        return None

    c_group = col("group_PDB")
    c_model = col("pdbx_PDB_model_num")
    c_chain = col("auth_asym_id", "label_asym_id")
    c_seq = col("auth_seq_id", "label_seq_id")
    c_icode = col("pdbx_PDB_ins_code")
    c_comp = col("auth_comp_id", "label_comp_id")
    c_name = col("auth_atom_id", "label_atom_id")
    c_alt = col("label_alt_id")
    c_occ = col("occupancy")
    c_b = col("B_iso_or_equiv")
    c_el = col("type_symbol")
    for needed, label in (
        (c_chain, "chaîne"),
        (c_seq, "numéro"),
        (c_comp, "résidu"),
        (c_name, "nom d'atome"),
    ):
        if needed is None:
            raise StructureFormatError(f"{origin} : colonne _atom_site de {label} absente")

    def value(row: gemmi.cif.Table.Row, index: int | None) -> str:
        if index is None:
            return ""
        raw = row.str(index)
        return "" if raw in ("?", ".") else raw

    atoms: list[RawAtom] = []
    for number, row in enumerate(table, start=1):
        where = f"{origin}, _atom_site ligne {number}"
        model_text = value(row, c_model)
        model = _strict_int(model_text, "numéro de modèle", where) if model_text else 1
        key: AtomKey = (
            model,
            value(row, c_chain),
            _strict_int(value(row, c_seq) or "?", "numéro de résidu", where),
            value(row, c_icode),
            value(row, c_comp),
            value(row, c_name),
            value(row, c_alt),
        )
        atoms.append(
            RawAtom(
                key,
                _strict_float(row.str(tags["Cartn_x"]), "coordonnée x", where),
                _strict_float(row.str(tags["Cartn_y"]), "coordonnée y", where),
                _strict_float(row.str(tags["Cartn_z"]), "coordonnée z", where),
                _optional_float(row.str(c_occ), "occupation", where) if c_occ is not None else None,
                _optional_float(row.str(c_b), "facteur B", where) if c_b is not None else None,
                value(row, c_el).upper() or None,
                value(row, c_group).upper() == "HETATM",
            )
        )
    return atoms


# --- Conversion depuis gemmi et comparaison à l'inventaire attendu -----------------------

_CONNECTION_KINDS = {
    gemmi.ConnectionType.Disulf: ConnectionKind.DISULFIDE,
    gemmi.ConnectionType.Covale: ConnectionKind.COVALENT,
    gemmi.ConnectionType.MetalC: ConnectionKind.METAL,
    gemmi.ConnectionType.Hydrog: ConnectionKind.HYDROGEN,
}

_COORDINATE_TOLERANCE = 1e-6
"""Écart admis entre la valeur relue par PepForge et celle de gemmi (arrondi flottant, Å)."""


def _residue_id(chain: str, residue: gemmi.Residue) -> ResidueId:
    number = residue.seqid.num
    if number is None:
        raise StructureFormatError(f"Résidu {residue.name} de la chaîne {chain} sans numéro")
    return ResidueId(chain, number, residue.seqid.icode.strip(), residue.name)


def _address_id(address: gemmi.AtomAddress) -> ResidueId:
    number = address.res_id.seqid.num
    if number is None:
        raise StructureFormatError("Liaison déclarée vers un résidu sans numéro")
    return ResidueId(
        address.chain_name, number, address.res_id.seqid.icode.strip(), address.res_id.name
    )


def _connections(structure: gemmi.Structure) -> tuple[Connection, ...]:
    return tuple(
        Connection(
            kind=_CONNECTION_KINDS.get(conn.type, ConnectionKind.UNKNOWN),
            first=_address_id(conn.partner1),
            first_atom=conn.partner1.atom_name,
            second=_address_id(conn.partner2),
            second_atom=conn.partner2.atom_name,
        )
        for conn in structure.connections
    )


def _key_of(model: int, chain: str, residue: gemmi.Residue, atom: gemmi.Atom) -> AtomKey:
    rid = _residue_id(chain, residue)
    altloc = "" if atom.altloc == "\0" else atom.altloc
    return (model, chain, rid.number, rid.insertion_code, residue.name, atom.name, altloc)


def _describe_key(key: AtomKey) -> str:
    model, chain, number, icode, residue, name, alt = key
    alt_text = f" conformation {alt}" if alt else ""
    return f"modèle {model}, {chain}:{residue}{number}{icode} atome {name}{alt_text}"


def _build_models(
    structure: gemmi.Structure, raw: list[RawAtom], origin: str
) -> tuple[dict[int, tuple[Atom, ...]], int, int]:
    """Atomes de gemmi, comparés un à un aux atomes relus dans le fichier brut."""
    expected: dict[AtomKey, RawAtom] = {}
    for item in raw:
        if item.key in expected:
            raise StructureFormatError(f"{origin} : atome en double ({_describe_key(item.key)})")
        expected[item.key] = item
    models: dict[int, tuple[Atom, ...]] = {}
    seen: set[AtomKey] = set()
    without_element = sum(1 for item in raw if item.element is None)
    without_occupancy = sum(1 for item in raw if item.occupancy is None)
    model_numbers = [model.num for model in structure]
    duplicates = [n for n, count in Counter(model_numbers).items() if count > 1]
    if duplicates:
        raise StructureFormatError(f"{origin} : numéros de modèle en double : {duplicates}")
    for model in structure:
        atoms: list[Atom] = []
        for chain in model:
            for residue in chain:
                rid = _residue_id(chain.name, residue)
                for atom in residue:
                    key = _key_of(model.num, chain.name, residue, atom)
                    source = expected.get(key)
                    if source is None or key in seen:
                        raise StructureFormatError(
                            f"{origin} : gemmi a produit un atome absent du fichier "
                            f"({_describe_key(key)})"
                        )
                    seen.add(key)
                    position = (atom.pos.x, atom.pos.y, atom.pos.z)
                    if any(
                        not math.isclose(a, b, abs_tol=_COORDINATE_TOLERANCE)
                        for a, b in zip(position, (source.x, source.y, source.z), strict=True)
                    ):
                        raise StructureFormatError(
                            f"{origin} : coordonnées de gemmi {position} différentes du fichier "
                            f"{(source.x, source.y, source.z)} ({_describe_key(key)})"
                        )
                    atoms.append(
                        Atom(
                            serial=atom.serial,
                            name=atom.name,
                            element=source.element or atom.element.name.upper(),
                            residue_name=residue.name,
                            residue_number=rid.number,
                            chain_id=chain.name,
                            x=source.x,
                            y=source.y,
                            z=source.z,
                            alt_loc=key[6] or None,
                            hetero=source.hetero,
                            insertion_code=rid.insertion_code,
                            occupancy=source.occupancy,
                            b_factor=source.b_factor,
                        )
                    )
        models[model.num] = tuple(atoms)
    lost = [key for key in expected if key not in seen]
    if lost:
        sample = "; ".join(_describe_key(k) for k in lost[:5])
        raise StructureFormatError(
            f"{origin} : {len(lost)} atome(s) du fichier absent(s) de la lecture par gemmi "
            f"({sample})"
        )
    return models, without_element, without_occupancy


# --- Résidus déclarés non observés --------------------------------------------------------

_REMARK_465 = re.compile(
    r"^\s*(?P<model>\d+)?\s*(?P<res>[A-Z0-9]{1,3})\s+(?P<chain>\S)\s+(?P<num>-?\d+)(?P<icode>[A-Z]?)\s*$"
)


def _remark_465(structure: gemmi.Structure, origin: str) -> dict[int | None, tuple[ResidueId, ...]]:
    found: dict[int | None, list[ResidueId]] = {}
    started = False
    for line in structure.raw_remarks:
        if not line.startswith("REMARK 465"):
            continue
        body = line[10:].rstrip()
        if not started:
            started = "RES C SSSEQI" in body
            continue
        if not body.strip():
            continue
        match = _REMARK_465.match(body)
        if not match:
            raise StructureFormatError(f"{origin} : ligne REMARK 465 illisible : {line.rstrip()!r}")
        model = int(match["model"]) if match["model"] else None
        found.setdefault(model, []).append(
            ResidueId(match["chain"], int(match["num"]), match["icode"], match["res"])
        )
    return {k: tuple(v) for k, v in found.items()}


def _unobserved_mmcif(
    block: gemmi.cif.Block, origin: str
) -> dict[int | None, tuple[ResidueId, ...]]:
    table = block.find_mmcif_category("_pdbx_unobs_or_zero_occ_residues.")
    if len(table) == 0:
        return {}
    tags = {
        tag.removeprefix("_pdbx_unobs_or_zero_occ_residues."): i for i, tag in enumerate(table.tags)
    }
    found: dict[int | None, list[ResidueId]] = {}
    for number, row in enumerate(table, start=1):
        where = f"{origin}, _pdbx_unobs_or_zero_occ_residues ligne {number}"

        def get(name: str, row: gemmi.cif.Table.Row = row) -> str:
            index = tags.get(name)
            if index is None:
                return ""
            raw = row.str(index)
            return "" if raw in ("?", ".") else raw

        # polymer_flag Y : résidu de polymère ; occupancy_flag 1 : non observé
        # (0 signalerait un résidu présent mais d'occupation nulle).
        if get("polymer_flag") != "Y" or get("occupancy_flag") != "1":
            continue
        model_text = get("PDB_model_num")
        model = _strict_int(model_text, "numéro de modèle", where) if model_text else None
        found.setdefault(model, []).append(
            ResidueId(
                get("auth_asym_id"),
                _strict_int(get("auth_seq_id") or "?", "numéro de résidu", where),
                get("PDB_ins_code"),
                get("auth_comp_id"),
            )
        )
    return {k: tuple(v) for k, v in found.items()}


# --- Lecteurs et écrivains inscrits au registre -------------------------------------------


def read_pdb(text: str, origin: str) -> ParsedFile:
    raw, declared = _scan_pdb(text, origin)
    if not raw:
        raise StructureFormatError(f"{origin} : aucune ligne ATOM ou HETATM")
    try:
        structure = gemmi.read_pdb_string(text)
    except (RuntimeError, ValueError) as exc:
        raise StructureFormatError(f"{origin} : lecture PDB refusée par gemmi ({exc})") from None
    models, no_element, no_occupancy = _build_models(structure, raw, origin)
    return ParsedFile(
        format=StructureFormat.PDB,
        models=models,
        connections=_connections(structure),
        declared_missing=_remark_465(structure, origin),
        declared_model_count=declared,
        atoms_without_element=no_element,
        atoms_without_occupancy=no_occupancy,
    )


def read_mmcif(text: str, origin: str) -> ParsedFile:
    try:
        document = gemmi.cif.read_string(text)
    except (RuntimeError, ValueError) as exc:
        raise StructureFormatError(f"{origin} : lecture mmCIF refusée par gemmi ({exc})") from None
    if len(document) != 1:
        raise StructureFormatError(
            f"{origin} : {len(document)} blocs de données ; "
            "un seul attendu (une structure par fichier)"
        )
    block = document.sole_block()
    raw = _scan_mmcif(block, origin)
    try:
        structure = gemmi.make_structure_from_block(block)
    except (RuntimeError, ValueError) as exc:
        raise StructureFormatError(f"{origin} : lecture mmCIF refusée par gemmi ({exc})") from None
    models, no_element, no_occupancy = _build_models(structure, raw, origin)
    return ParsedFile(
        format=StructureFormat.MMCIF,
        models=models,
        connections=_connections(structure),
        declared_missing=_unobserved_mmcif(block, origin),
        declared_model_count=None,
        atoms_without_element=no_element,
        atoms_without_occupancy=no_occupancy,
    )


def _to_gemmi(atoms: tuple[Atom, ...], connections: tuple[Connection, ...]) -> gemmi.Structure:
    """Une seule structure gemmi, un modèle, chaînes et résidus dans l'ordre d'apparition."""
    structure = gemmi.Structure()
    model = gemmi.Model(1)
    chains: dict[str, dict[ResidueId, list[Atom]]] = {}
    for entry in atoms:
        chains.setdefault(entry.chain_id, {}).setdefault(entry.residue_id, []).append(entry)
    for chain_id, residues in chains.items():
        chain = gemmi.Chain(chain_id)
        for rid, members in residues.items():
            residue = gemmi.Residue()
            residue.name = rid.name
            residue.seqid = gemmi.SeqId(rid.number, rid.insertion_code or " ")
            residue.het_flag = "H" if members[0].hetero else "A"
            for item in members:
                g_atom = gemmi.Atom()
                g_atom.name = item.name
                g_atom.element = gemmi.Element(item.element)
                g_atom.pos = gemmi.Position(item.x, item.y, item.z)
                g_atom.altloc = item.alt_loc or "\0"
                g_atom.occ = item.occupancy if item.occupancy is not None else 1.0
                g_atom.b_iso = item.b_factor if item.b_factor is not None else 0.0
                g_atom.serial = item.serial
                residue.add_atom(g_atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    for index, link in enumerate(connections, start=1):
        conn = gemmi.Connection()
        conn.name = f"conn{index}"
        conn.type = next(
            (t for t, k in _CONNECTION_KINDS.items() if k == link.kind),
            gemmi.ConnectionType.Unknown,
        )
        conn.partner1 = gemmi.AtomAddress(
            link.first.chain_id,
            gemmi.SeqId(link.first.number, link.first.insertion_code or " "),
            link.first.name,
            link.first_atom,
        )
        conn.partner2 = gemmi.AtomAddress(
            link.second.chain_id,
            gemmi.SeqId(link.second.number, link.second.insertion_code or " "),
            link.second.name,
            link.second_atom,
        )
        structure.connections.append(conn)
    structure.setup_cell_images()
    return structure


def write_pdb_text(atoms: tuple[Atom, ...], connections: tuple[Connection, ...]) -> str:
    return _to_gemmi(atoms, connections).make_pdb_string(gemmi.PdbWriteOptions())


def write_mmcif_text(atoms: tuple[Atom, ...], connections: tuple[Connection, ...]) -> str:
    structure = _to_gemmi(atoms, connections)
    structure.name = "pepforge"
    return structure.make_mmcif_document().as_string()


Reader = Callable[[str, str], ParsedFile]
Writer = Callable[[tuple[Atom, ...], tuple[Connection, ...]], str]

READERS: dict[StructureFormat, Reader] = {
    StructureFormat.PDB: read_pdb,
    StructureFormat.MMCIF: read_mmcif,
}
WRITERS: dict[StructureFormat, Writer] = {
    StructureFormat.PDB: write_pdb_text,
    StructureFormat.MMCIF: write_mmcif_text,
}
