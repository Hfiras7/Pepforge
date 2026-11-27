"""Lecture et écriture de structures 3D, depuis un fichier local ou une URL.

Deux modes d'importation, une seule lecture (demande de Firas, 27 novembre 2026) :

- `read_structure(chemin, ...)` lit un fichier du disque ;
- `read_structure_url(url, ...)` le télécharge d'abord sous contrôle
  (`pepforge.io.retrieval.download_file`), puis le lit de la même façon.

Ni l'un ni l'autre ne corrige quoi que ce soit (architecture, section 5.7) :

- **format** : annoncé par l'extension, reconnu sur le contenu, ou imposé par
  `fmt` ; s'ils se contredisent, `StructureFormatError`. Un contenu qui ne
  ressemble à aucun format pris en charge est refusé ;
- **modèles** : un fichier multi-modèles (RMN) exige `model=` ; le premier
  n'est jamais pris en silence. Un en-tête `NUMMDL` qui annonce un autre nombre
  de modèles que le fichier n'en contient signale un fichier tronqué : refus ;
- **conformations alternatives, codes d'insertion, hydrogènes, eaux, ligands** :
  tout est gardé tel quel. Leur résolution est une décision, prise plus tard par
  un `StructureModifier` journalisé (`pepforge.build`) ;
- **liaisons déclarées** (`SSBOND`, `LINK`, `struct_conn`) et **résidus déclarés
  non observés** (`REMARK 465`, `pdbx_unobs_or_zero_occ_residues`) : gardés
  dans `Structure`, pour que rien de ce que le fichier annonce ne disparaisse.

Chaque lecture est inscrite au journal avec sa provenance (chemin ou URL, date,
taille, SHA-256), le format, le modèle retenu et la version de gemmi.

L'écriture (`write_structure`) relit le fichier qu'elle vient d'écrire et le
compare à la structure, atome par atome ; le fichier n'apparaît à son chemin
final qu'une fois cette comparaison réussie.
"""

from __future__ import annotations

import math
import os
import tempfile
from collections import Counter
from collections.abc import Mapping
from pathlib import Path

from pepforge.core.structure import (
    STRUCTURE_DECLARED_MODEL_COUNT,
    STRUCTURE_FORMAT,
    STRUCTURE_MODEL_COUNT,
    STRUCTURE_SHA256,
    Atom,
    ConnectionKind,
    Structure,
)
from pepforge.exceptions import StructureFormatError, UnsupportedFeatureError
from pepforge.io.retrieval import DownloadLimits, RetrievedFile, download_file, load_local_file
from pepforge.io.structure_formats import (
    GEMMI_VERSION,
    READERS,
    WRITERS,
    ParsedFile,
    StructureFormat,
    decode_text,
    detect_format,
    format_from_name,
)
from pepforge.report.journal import Journal

PDB_MAX_SERIAL = 99_999
PDB_CONNECTION_KINDS = frozenset(
    {ConnectionKind.DISULFIDE, ConnectionKind.COVALENT, ConnectionKind.METAL}
)
PDB_RESIDUE_NUMBER_RANGE = (-999, 9_999)

_WRITE_TOLERANCE: Mapping[StructureFormat, float] = {
    # Le PDB écrit les coordonnées avec trois décimales : écart d'arrondi ≤ 0,0005 Å.
    StructureFormat.PDB: 5e-4 + 1e-9,
    StructureFormat.MMCIF: 5e-4 + 1e-9,
}
_OCCUPANCY_TOLERANCE = 5e-3 + 1e-9
"""Le PDB écrit occupation et facteur B avec deux décimales."""


def _as_format(value: StructureFormat | str | None) -> StructureFormat | None:
    if value is None:
        return None
    try:
        return StructureFormat(value)
    except ValueError:
        known = ", ".join(f.value for f in StructureFormat)
        raise StructureFormatError(
            f"Format inconnu : {value!r} (formats pris en charge : {known})"
        ) from None


def resolve_format(
    text: str, name: str, origin: str, declared: StructureFormat | str | None = None
) -> StructureFormat:
    """Format du fichier : extension, contenu et format imposé doivent concorder."""
    requested = _as_format(declared)
    by_name = format_from_name(name)
    by_content = detect_format(text)
    if by_content is None:
        raise StructureFormatError(
            f"{origin} : contenu non reconnu comme {', '.join(f.value for f in StructureFormat)}"
        )
    for label, value in (("l'extension", by_name), ("le format demandé", requested)):
        if value is not None and value != by_content:
            raise StructureFormatError(
                f"{origin} : {label} annonce {value.value}, "
                f"mais le contenu est du {by_content.value}"
            )
    return by_content


def _select_model(parsed: ParsedFile, model: int | None, origin: str) -> int:
    numbers = sorted(parsed.models)
    if parsed.declared_model_count is not None and parsed.declared_model_count != len(numbers):
        raise StructureFormatError(
            f"{origin} : l'en-tête NUMMDL annonce {parsed.declared_model_count} modèle(s), "
            f"le fichier en contient {len(numbers)} (fichier tronqué ou en-tête incohérent)"
        )
    if model is None:
        if len(numbers) > 1:
            raise StructureFormatError(
                f"{origin} : {len(numbers)} modèles ({numbers[0]} à {numbers[-1]}) ; "
                "choisir le modèle à lire (model= dans l'API, --model en ligne de commande) ; "
                "aucun n'est pris par défaut"
            )
        return numbers[0]
    if model not in parsed.models:
        raise StructureFormatError(
            f"{origin} : modèle {model} absent ; modèles présents : {numbers}"
        )
    return model


def parse_structure(
    retrieved: RetrievedFile,
    *,
    journal: Journal,
    model: int | None = None,
    fmt: StructureFormat | str | None = None,
) -> Structure:
    """Lit un fichier déjà obtenu (local ou téléchargé) et journalise la lecture."""
    origin = retrieved.origin
    text = decode_text(retrieved.content, origin)
    structure_format = resolve_format(text, retrieved.name, origin, fmt)
    by_extension = format_from_name(retrieved.name)
    parsed = READERS[structure_format](text, origin)
    chosen = _select_model(parsed, model, origin)
    atoms = parsed.models[chosen]
    if not atoms:
        raise StructureFormatError(f"{origin} : le modèle {chosen} ne contient aucun atome")
    missing = parsed.declared_missing.get(None, ()) + parsed.declared_missing.get(chosen, ())
    missing_atoms = parsed.declared_missing_atoms.get(None, ()) + parsed.declared_missing_atoms.get(
        chosen, ()
    )
    metadata = {
        STRUCTURE_FORMAT: structure_format.value,
        STRUCTURE_SHA256: retrieved.sha256,
        STRUCTURE_MODEL_COUNT: str(len(parsed.models)),
    }
    if parsed.declared_model_count is not None:
        metadata[STRUCTURE_DECLARED_MODEL_COUNT] = str(parsed.declared_model_count)
    structure = Structure(
        atoms=atoms,
        source=origin,
        metadata=metadata,
        model=chosen if len(parsed.models) > 1 or model is not None else None,
        connections=parsed.connections,
        declared_missing_residues=missing,
        declared_missing_atoms=missing_atoms,
    )
    journal.record(
        "read_structure",
        "importation d'une structure, sans aucune correction",
        **retrieved.provenance(),
        format=structure_format.value,
        format_from_extension=by_extension.value if by_extension else None,
        model=chosen,
        model_count=len(parsed.models),
        atoms=len(atoms),
        atoms_without_element=parsed.atoms_without_element,
        atoms_without_occupancy=parsed.atoms_without_occupancy,
        declared_connections=len(parsed.connections),
        declared_missing_residues=len(missing),
        declared_missing_atoms=len(missing_atoms),
        reader=f"gemmi {GEMMI_VERSION}",
    )
    return structure


def read_structure(
    path: str | Path,
    *,
    journal: Journal,
    model: int | None = None,
    fmt: StructureFormat | str | None = None,
    max_bytes: int = DownloadLimits.max_bytes,
) -> Structure:
    """Importe une structure depuis un fichier local (PDB, mmCIF, éventuellement compressé gzip)."""
    return parse_structure(
        load_local_file(path, max_bytes=max_bytes), journal=journal, model=model, fmt=fmt
    )


def read_structure_url(
    url: str,
    *,
    journal: Journal,
    model: int | None = None,
    fmt: StructureFormat | str | None = None,
    limits: DownloadLimits | None = None,
) -> Structure:
    """Importe une structure depuis une URL, par un téléchargement contrôlé puis la même lecture."""
    return parse_structure(download_file(url, limits), journal=journal, model=model, fmt=fmt)


# --- Écriture ------------------------------------------------------------------------------


def _pdb_limits(structure: Structure) -> list[str]:
    problems: list[str] = []
    if len(structure.atoms) > PDB_MAX_SERIAL:
        problems.append(f"{len(structure.atoms)} atomes, le PDB en admet {PDB_MAX_SERIAL}")
    long_chains = sorted({a.chain_id for a in structure.atoms if len(a.chain_id) > 1})
    if long_chains:
        problems.append(f"identifiants de chaîne de plus d'un caractère : {long_chains}")
    low, high = PDB_RESIDUE_NUMBER_RANGE
    out = sorted({a.residue_number for a in structure.atoms if not low <= a.residue_number <= high})
    if out:
        problems.append(f"numéros de résidu hors de [{low}, {high}] : {out[:5]}")
    long_names = sorted({a.residue_name for a in structure.atoms if len(a.residue_name) > 3})
    if long_names:
        problems.append(f"noms de résidu de plus de trois caractères : {long_names}")
    long_atoms = sorted({a.name for a in structure.atoms if len(a.name) > 4})
    if long_atoms:
        problems.append(f"noms d'atome de plus de quatre caractères : {long_atoms}")
    # Le PDB porte les ponts disulfure (SSBOND) et les liaisons covalentes ou de
    # coordination (LINK) ; il n'a pas d'enregistrement pour une liaison hydrogène.
    unwritable = Counter(
        c.kind for c in structure.connections if c.kind not in PDB_CONNECTION_KINDS
    )
    if unwritable:
        detail = ", ".join(f"{k} ×{v}" for k, v in sorted(unwritable.items()))
        problems.append(f"liaisons déclarées sans équivalent PDB : {detail}")
    return problems


def _compare_written(
    expected: tuple[Atom, ...], written: tuple[Atom, ...], fmt: StructureFormat, where: str
) -> None:
    """Inventaire attendu (la structure) contre fichier produit (relu), atome par atome."""

    def key(atom: Atom) -> tuple[str, int, str, str, str, str]:
        return (
            atom.chain_id,
            atom.residue_number,
            atom.insertion_code,
            atom.residue_name,
            atom.name,
            atom.alt_loc or "",
        )

    produced = {key(a): a for a in written}
    if len(produced) != len(written) or len(written) != len(expected):
        raise StructureFormatError(
            f"{where} : {len(written)} atomes relus pour {len(expected)} écrits"
        )
    tolerance = _WRITE_TOLERANCE[fmt]
    for atom in expected:
        other = produced.get(key(atom))
        if other is None:
            raise StructureFormatError(f"{where} : atome {key(atom)} absent du fichier écrit")
        position = (atom.x, atom.y, atom.z)
        read_back = (other.x, other.y, other.z)
        if any(abs(a - b) > tolerance for a, b in zip(position, read_back, strict=True)):
            raise StructureFormatError(
                f"{where} : atome {key(atom)} relu en {read_back} au lieu de {position}"
            )
        if other.element != atom.element or other.hetero != atom.hetero:
            raise StructureFormatError(
                f"{where} : élément ou type d'enregistrement de {key(atom)} perdu"
            )
        for label, original, back in (
            ("occupation", atom.occupancy, other.occupancy),
            ("facteur B", atom.b_factor, other.b_factor),
        ):
            if original is not None and (
                back is None or not math.isclose(original, back, abs_tol=_OCCUPANCY_TOLERANCE)
            ):
                raise StructureFormatError(
                    f"{where} : {label} de {key(atom)} relue {back}, {original} écrite"
                )


def write_structure(
    structure: Structure,
    path: str | Path,
    fmt: StructureFormat | str,
    *,
    journal: Journal,
) -> Path:
    """Écrit `structure` au format `fmt`, relit le fichier et le compare avant de le livrer.

    Une structure que le format ne peut pas représenter (au-delà de 99 999
    atomes ou chaîne de deux caractères en PDB) lève `UnsupportedFeatureError`
    au lieu d'être tronquée. Une occupation ou un facteur B absents sont écrits
    1,00 et 0,00, comme le format l'exige ; leur nombre est inscrit au journal.
    """
    target = Path(path)
    structure_format = _as_format(fmt)
    assert structure_format is not None
    if not structure.atoms:
        raise StructureFormatError("Structure vide : rien à écrire")
    if structure_format is StructureFormat.PDB:
        problems = _pdb_limits(structure)
        if problems:
            raise UnsupportedFeatureError(
                "Structure non représentable en PDB sans perte : "
                + "; ".join(problems)
                + " ; écrire en mmCIF"
            )
    text = WRITERS[structure_format](structure.atoms, structure.connections)
    target.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(
        dir=target.parent, prefix=f".{target.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
        # Le fichier réellement écrit sur le disque est relu, pas le texte en mémoire.
        written_text = decode_text(Path(temporary).read_bytes(), str(target))
        reread = READERS[structure_format](written_text, str(target))
        if len(reread.models) != 1:
            raise StructureFormatError(f"{target} : {len(reread.models)} modèles écrits, 1 attendu")
        (written,) = reread.models.values()
        _compare_written(structure.atoms, written, structure_format, str(target))
        if len(reread.connections) != len(structure.connections):
            raise StructureFormatError(
                f"{target} : {len(reread.connections)} liaisons déclarées relues pour "
                f"{len(structure.connections)} écrites"
            )
        os.replace(temporary, target)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    journal.record(
        "write_structure",
        "écriture d'une structure, relue et comparée atome par atome",
        path=str(target.resolve()),
        format=structure_format.value,
        atoms=len(structure.atoms),
        connections=len(structure.connections),
        occupancy_written_as_default=sum(1 for a in structure.atoms if a.occupancy is None),
        b_factor_written_as_default=sum(1 for a in structure.atoms if a.b_factor is None),
        writer=f"gemmi {GEMMI_VERSION}",
    )
    return target
