"""Validation du module 3 sur des structures réelles de la PDB.

Ce script n'est pas un test pytest : il demande le réseau. Il se relance tel quel :

    pip install -e ".[validation]"
    python validation/module03_structures.py

Source des structures. Le réseau de l'environnement de développement refusait
les serveurs de la PDB (files.rcsb.org, www.ebi.ac.uk : HTTP 403 du proxy, 27
novembre 2026). Les fichiers sont donc pris dans la suite de tests de
Biopython 1.88, qui distribue des copies d'entrées de la PDB, à l'étiquette
Git `biopython-188` (adresse figée). Le fichier `1A8O.cif` de cette étiquette a
été vérifié identique, octet pour octet, à celui de l'archive source publiée
sur PyPI (archive biopython-1.88.tar.gz, SHA-256 9aaa31c0...071a, conforme à
l'empreinte affichée par PyPI). Certaines copies ont été raccourcies par
Biopython (un seul modèle gardé sur dix) : c'est voulu ici, c'est le cas du
fichier tronqué que le lecteur doit refuser.

Chaque fichier est importé **par URL** avec `read_structure_url` (téléchargement
contrôlé réel, en https), puis comparé à une lecture indépendante :

1. **Biopython 1.88** (Hamelryck et Manderick, Bioinformatics 19:2308, 2003),
   autre lecteur publié, en mode strict pour le PDB (`PERMISSIVE=0`) : nombre de
   modèles, et pour **chaque** modèle l'ensemble des atomes (chaîne, numéro,
   code d'insertion, résidu, nom, conformation, coordonnées à 0,001 Å) ;
2. **l'en-tête du fichier lui-même**, lu par `MMCIF2Dict` ou ligne à ligne pour
   le PDB : résidus déclarés non observés et liaisons déclarées ;
3. **le même fichier dans l'autre format** quand la PDB fournit les deux ;
4. **l'écriture puis la relecture** en PDB et en mmCIF.

Les résultats sont écrits dans `validation/results/module03_structures.yaml`.
"""

from __future__ import annotations

import sys
import warnings
from collections.abc import Iterable
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import Bio
import yaml
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from pepforge.core.peptide import C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED
from pepforge.core.structure import Atom
from pepforge.exceptions import ExtractionRefusedError, PepForgeError, StructureError
from pepforge.io import (
    ExtractionReport,
    extract_peptides,
    inventory,
    read_structure,
    read_structure_url,
)
from pepforge.io.retrieval import download_file
from pepforge.io.structure import write_structure
from pepforge.io.structure_formats import GEMMI_VERSION
from pepforge.report.journal import Journal

BASE = "https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/PDB/"

# Fichier → ce qu'il exerce. Aucun n'est une cible de la thèse : ce sont des cas de lecture.
FILES: dict[str, str] = {
    "1A8O.pdb": "rayons X ; sélénométhionine MSE en HETATM ; eaux ; pont disulfure ; LINK",
    "1A8O.cif": "même entrée en mmCIF : MSE y est en ATOM",
    "2XHE.pdb": "complexe de deux protéines (6 315 atomes) ; 142 résidus non observés",
    "2XHE.cif": "même entrée en mmCIF",
    "1LCD.pdb": "RMN, 3 modèles ; complexe protéine–ADN ; ion Na",
    "1LCD.cif": "même entrée en mmCIF ; liaisons hydrogène déclarées dans struct_conn",
    "2BEG.cif": "RMN, 10 modèles ; fibrille Aβ(17-42), 5 chaînes peptidiques ; hydrogènes",
    "2BEG.pdb": "copie raccourcie par Biopython : 1 modèle, mais NUMMDL 10",
    "1AS5.cif": "RMN, 14 modèles ; petit peptide à hydroxyprolines HYP et C-ter NH2",
    "2n0n_M1.pdb": "copie d'un modèle d'un ensemble RMN : NUMMDL 20 ; AIB, NH2",
    "7DDO.pdb": "cryo-ME ; conformations alternatives ; ion Zn ; glycanes NAG",
    "6WQA.cif": "98 atomes à conformations alternatives ; ligand",
    "4ZHL.cif": "codes d'insertion (numérotation chymotrypsine) ; chaîne peptidique P de 10 résidus",
    "4CUP.cif": "conformations alternatives ; ligands ; résidus non observés",
}


# Motifs de refus attendus pour chaque chaîne peptidique, établis à la main dans l'analyse
# du 27 novembre 2026 (bilans/module-03-extraction-analyse.md, section 1) avant le code :
# liste vide = extraction acceptée. Les chaînes d'ADN de 1LCD (B, C) sont refusées comme
# non peptidiques. 2BEG.pdb et 2n0n_M1.pdb ne sont pas lus (NUMMDL incohérent), donc absents.
EXPECTED_REFUSALS: dict[tuple[str, str], list[str]] = {
    **{(f, "A"): ["modified_residue"] for f in ("1A8O.pdb", "1A8O.cif", "1AS5.cif")},
    **{(f, "A"): ["break_declared", "unobserved_C"] for f in ("2XHE.pdb", "2XHE.cif")},
    **{
        (f, "B"): ["break_declared", "unobserved_C", "unobserved_N"]
        for f in ("2XHE.pdb", "2XHE.cif")
    },
    **{(f, "A"): [] for f in ("1LCD.pdb", "1LCD.cif")},
    **{("2BEG.cif", c): ["unobserved_N"] for c in "ABCDE"},
    ("7DDO.pdb", "A"): ["side_chain_modification"],
    ("7DDO.pdb", "C"): ["side_chain_modification", "unobserved_C", "unobserved_N"],
    ("6WQA.cif", "A"): ["break_undeclared"],
    ("4ZHL.cif", "U"): [],
    ("4ZHL.cif", "P"): [],
    ("4CUP.cif", "A"): ["unobserved_C"],
}


def _refusal_causes(report: ExtractionReport) -> list[str]:
    """Motifs de refus d'un rapport, par catégorie."""
    causes = {j.category.value for j in report.junctions if j.decision.value == "reject"}
    causes |= {
        f"unobserved_{t.side}" for t in report.unobserved_termini if t.decision.value == "reject"
    }
    causes |= {c.category.value for c in report.components if c.decision.value == "reject"}
    for word, cause in (
        ("sans place", "unplaced_missing"),
        ("liaison déclarée", "declared_link"),
        ("disulfure", "disulfide"),
        ("microhétérogénéité", "microheterogeneity"),
        ("groupe terminal", "terminal_cap"),
    ):
        if any(word in reason for reason in report.refusals):
            causes.add(cause)
    return sorted(causes)


def _report_summary(report: ExtractionReport) -> dict[str, Any]:
    return {
        "accepted": report.accepted,
        "junctions_by_category": report.junctions_by_category(),
        "notable_junctions": [
            f"{j.before} → {j.after} : {j.category.value}, "
            + ", ".join(f"{m.distance:.3f} Å" for m in j.measurements)
            for j in report.junctions
            if j.category.value not in ("continuous",) and not j.insertion_code
        ],
        "insertion_code_junctions": sum(1 for j in report.junctions if j.insertion_code),
        "components": {
            c.residue: c.category.value for c in report.components if c.category.value != "water"
        },
        "waters": sum(1 for c in report.components if c.category.value == "water"),
        "unobserved_termini": {t.side: len(t.residues) for t in report.unobserved_termini},
        "declared_missing_atoms": len(report.declared_missing_atoms),
        "absent_backbone_atoms": list(report.absent_backbone_atoms),
        "refusal_causes": _refusal_causes(report),
    }


def _pepforge_key(atom: Atom) -> tuple[Any, ...]:
    return (
        atom.chain_id,
        atom.residue_number,
        atom.insertion_code,
        atom.residue_name,
        atom.name,
        atom.alt_loc or "",
        round(atom.x, 3),
        round(atom.y, 3),
        round(atom.z, 3),
    )


def _biopython_atoms(model: Any) -> set[tuple[Any, ...]]:
    keys = set()
    for chain in model:
        for residue in chain:
            _, number, icode = residue.id
            for atom in residue.get_unpacked_list():
                x, y, z = (round(float(v), 3) for v in atom.coord)
                keys.add(
                    (
                        chain.id,
                        number,
                        icode.strip(),
                        residue.get_resname(),
                        atom.get_name(),
                        atom.get_altloc().strip(),
                        x,
                        y,
                        z,
                    )
                )
    return keys


def _biopython(path: Path) -> tuple[Any, list[str]]:
    caught: list[str] = []
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        if path.suffix == ".cif":
            structure = MMCIFParser(QUIET=False).get_structure("x", str(path))
        else:
            structure = PDBParser(PERMISSIVE=0, QUIET=False).get_structure("x", str(path))
    caught = sorted({str(w.message)[:120] for w in record})
    return structure, caught


def _declared_from_header(path: Path) -> dict[str, int]:
    """Résidus non observés et liaisons déclarés, comptés sans PepForge ni gemmi."""
    if path.suffix == ".cif":
        data = MMCIF2Dict(str(path))
        flags = zip(
            data.get("_pdbx_unobs_or_zero_occ_residues.polymer_flag", []),
            data.get("_pdbx_unobs_or_zero_occ_residues.occupancy_flag", []),
            data.get("_pdbx_unobs_or_zero_occ_residues.PDB_model_num", []),
            strict=True,
        )
        unobserved_model_1 = sum(1 for p, o, m in flags if p == "Y" and o == "1" and m == "1")
        connections = len(data.get("_struct_conn.id", []))
        return {"unobserved_model_1": unobserved_model_1, "connections": connections}
    lines = path.read_text().splitlines()
    remark = [line for line in lines if line.startswith("REMARK 465")]
    start = next((i for i, line in enumerate(remark) if "RES C SSSEQI" in line), None)
    unobserved = 0 if start is None else sum(1 for line in remark[start + 1 :] if line[10:].strip())
    connections = sum(1 for line in lines if line.startswith(("SSBOND", "LINK  ")))
    return {"unobserved_model_1": unobserved, "connections": connections}


def _entity_sequence(path: Path, chain: str) -> str | None:
    """Séquence canonique déclarée par le fichier mmCIF pour la chaîne (auteur) `chain`."""
    if path.suffix != ".cif":
        return None
    data = MMCIF2Dict(str(path))
    strands = data.get("_entity_poly.pdbx_strand_id", [])
    codes = data.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
    for strand, code in zip(strands, codes, strict=True):
        if chain in strand.split(","):
            return "".join(code.split())
    return None


def _models_of(structure_path: Path) -> list[int]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bio, _ = _biopython(structure_path)
    return [model.serial_num if hasattr(model, "serial_num") else model.id + 1 for model in bio]


def _check(results: dict[str, Any], name: str, ok: bool, detail: str) -> None:
    results.setdefault("checks", []).append({"check": name, "ok": bool(ok), "detail": detail})


def validate(name: str, workdir: Path) -> dict[str, Any]:
    url = BASE + name
    result: dict[str, Any] = {"file": name, "url": url, "exercises": FILES[name]}
    retrieved = download_file(url)
    local = workdir / name
    local.write_bytes(retrieved.content)
    result["sha256"] = retrieved.sha256
    result["size_bytes"] = retrieved.size

    bio, bio_warnings = _biopython(local)
    bio_models = list(bio)
    result["biopython_warnings"] = bio_warnings
    numbers = [m.serial_num for m in bio_models]

    # Lecture par URL : modèle 1 explicite si l'ensemble en a plusieurs.
    journal = Journal()
    try:
        first = read_structure_url(
            url, journal=journal, model=numbers[0] if len(numbers) > 1 else None
        )
    except StructureError as exc:
        result["pepforge_read"] = f"REFUSÉ : {exc}"
        _check(
            result,
            "refus attendu d'un fichier tronqué",
            "NUMMDL" in str(exc),
            f"Biopython lit le fichier sans erreur ({len(numbers)} modèle(s)) ; PepForge refuse.",
        )
        return result
    params = journal.entries[0].parameters
    _check(
        result,
        "provenance au journal",
        params["sha256"] == retrieved.sha256 and params["mode"] == "url",
        f"sha256 {params['sha256'][:12]}..., {params['size_bytes']} octets, {params['reader']}",
    )

    # 1. Biopython, modèle par modèle.
    _check(
        result,
        "nombre de modèles = Biopython",
        int(first.metadata["structure_model_count"]) == len(bio_models),
        f"PepForge {first.metadata['structure_model_count']}, Biopython {len(bio_models)}",
    )
    mismatches = []
    total = 0
    for bio_model in bio_models:
        model_number = bio_model.serial_num
        ours = read_structure(
            local, journal=Journal(), model=model_number if len(numbers) > 1 else None
        )
        mine = {_pepforge_key(a) for a in ours.atoms}
        theirs = _biopython_atoms(bio_model)
        total += len(ours.atoms)
        if mine != theirs or len(mine) != len(ours.atoms):
            mismatches.append(
                f"modèle {model_number} : {len(mine - theirs)} en trop, {len(theirs - mine)} manquants"
            )
    _check(
        result,
        "atomes identiques à Biopython, pour chaque modèle",
        not mismatches,
        "; ".join(mismatches) or f"{total} atomes comparés sur {len(bio_models)} modèle(s)",
    )

    # 2. En-tête du fichier.
    header = _declared_from_header(local)
    declared_missing = [r for r in first.declared_missing_residues]
    _check(
        result,
        "résidus non observés = en-tête",
        len(declared_missing) == header["unobserved_model_1"],
        f"PepForge {len(declared_missing)}, en-tête {header['unobserved_model_1']}",
    )
    _check(
        result,
        "liaisons déclarées = en-tête",
        len(first.connections) == header["connections"],
        f"PepForge {len(first.connections)}, en-tête {header['connections']}",
    )

    # 3. Inventaire (consigné pour relecture).
    inv = inventory(first)
    result["inventory"] = {
        "atoms": inv.atoms,
        "hydrogens": inv.hydrogens,
        "chains": {
            c.chain_id: {
                "standard": c.standard_residues,
                "amino_acid_like": c.amino_acid_like_residues,
                "waters": c.waters,
                "ions": list(c.ions),
                "unrecognized": c.unrecognized,
                "numbering_gaps": [list(g) for g in c.numbering_gaps],
                "insertion_codes": len(c.insertion_codes),
            }
            for c in inv.chains
        },
        "alternate_location_residues": len(inv.alternate_locations),
        "declared_missing_residues": len(inv.declared_missing_residues),
        "connections": inv.connections,
        "pending_decisions": list(inv.pending_decisions),
    }
    bio_waters = sum(1 for r in bio_models[0].get_residues() if r.id[0] == "W")
    _check(
        result,
        "eaux = Biopython",
        sum(c.waters for c in inv.chains) == bio_waters,
        f"PepForge {sum(c.waters for c in inv.chains)}, Biopython {bio_waters}",
    )
    bio_alt = sum(
        1
        for r in bio_models[0].get_residues()
        if any(a.get_altloc().strip() for a in r.get_unpacked_list())
    )
    _check(
        result,
        "résidus à conformations alternatives = Biopython",
        len(inv.alternate_locations) == bio_alt,
        f"PepForge {len(inv.alternate_locations)}, Biopython {bio_alt}",
    )

    # 4. Écriture et relecture.
    for fmt in ("pdb", "mmcif"):
        target = workdir / f"roundtrip_{name}.{fmt}"
        try:
            write_structure(first, target, fmt, journal=Journal())
            again = read_structure(target, journal=Journal())
            same = sorted(map(_pepforge_key, again.atoms)) == sorted(
                map(_pepforge_key, first.atoms)
            )
            _check(
                result, f"écriture {fmt} puis relecture sans perte", same, f"{len(again)} atomes"
            )
        except PepForgeError as exc:
            _check(
                result,
                f"écriture {fmt} refusée plutôt qu'incomplète",
                "sans équivalent PDB" in str(exc),
                str(exc)[:200],
            )

    # 5. Extraction d'un peptide, chaîne par chaîne : rapport et motifs de refus attendus.
    peptides = {}
    for chain in first.chain_ids:
        try:
            extraction = extract_peptides(first, chain, journal=Journal())
        except ExtractionRefusedError as exc:
            assert isinstance(exc.report, ExtractionReport)
            report = exc.report
            peptides[chain] = _report_summary(report)
        except PepForgeError as exc:
            peptides[chain] = f"refusé ({type(exc).__name__}) : {str(exc)[:160]}"
            _check(
                result,
                f"chaîne {chain} refusée comme non peptidique",
                (name, chain) not in EXPECTED_REFUSALS and "non peptidique" in str(exc),
                str(exc)[:160],
            )
            continue
        else:
            report = extraction.report
            (peptide,) = extraction.peptides
            peptides[chain] = {
                **_report_summary(report),
                "residues": "-".join(r.code for r in peptide.residues),
                "n_terminus_from_structure": N_TERMINUS_DEFAULTED not in peptide.metadata,
                "c_terminus_from_structure": C_TERMINUS_DEFAULTED not in peptide.metadata,
                "n_terminus": peptide.n_terminus.key,
                "c_terminus": peptide.c_terminus.key,
                "disulfides": [list(link.positions) for link in peptide.links],
            }
            declared = _entity_sequence(local, chain)
            if declared is not None:
                read = "".join(r.one_letter or "?" for r in peptide.residues)
                _check(
                    result,
                    f"séquence de la chaîne {chain} = séquence déclarée (_entity_poly)",
                    read == declared,
                    f"lue {read}, déclarée {declared}",
                )
        expected = EXPECTED_REFUSALS.get((name, chain))
        causes = _refusal_causes(report)
        _check(
            result,
            f"chaîne {chain} : motifs de refus attendus par l'analyse",
            expected is not None and causes == expected,
            f"obtenus {causes or 'aucun (acceptée)'}, attendus {expected}",
        )
    result["peptides"] = peptides
    return result


def _cross_format(results: dict[str, dict[str, Any]], workdir: Path) -> list[dict[str, Any]]:
    out = []
    for entry in ("1A8O", "2XHE", "1LCD"):
        pdb = workdir / f"{entry}.pdb"
        cif = workdir / f"{entry}.cif"
        model = 1 if entry == "1LCD" else None
        a = read_structure(pdb, journal=Journal(), model=model)
        b = read_structure(cif, journal=Journal(), model=model)
        same_atoms = sorted(map(_pepforge_key, a.atoms)) == sorted(map(_pepforge_key, b.atoms))
        out.append(
            {
                "entry": entry,
                "same_atoms_pdb_vs_mmcif": same_atoms,
                "atoms": [len(a), len(b)],
                "connections_pdb_vs_mmcif": [len(a.connections), len(b.connections)],
                "hetero_records_pdb_vs_mmcif": [
                    sum(x.hetero for x in a.atoms),
                    sum(x.hetero for x in b.atoms),
                ],
            }
        )
    return out


def main(names: Iterable[str] = FILES) -> int:
    report: dict[str, Any] = {
        "pepforge_reader": f"gemmi {GEMMI_VERSION}",
        "reference_reader": f"Biopython {Bio.__version__}",
        "source": BASE,
        "files": [],
    }
    with TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        results = {}
        for name in names:
            print(f"… {name}", file=sys.stderr)
            results[name] = validate(name, workdir)
            report["files"].append(results[name])
        report["cross_format"] = _cross_format(results, workdir)
    failed = [
        (f["file"], c["check"]) for f in report["files"] for c in f.get("checks", []) if not c["ok"]
    ]
    failed += [
        (c["entry"], "pdb = mmcif")
        for c in report["cross_format"]
        if not c["same_atoms_pdb_vs_mmcif"]
    ]
    report["summary"] = {
        "checks": sum(len(f.get("checks", [])) for f in report["files"]),
        "failed": [list(f) for f in failed],
    }
    out = Path(__file__).parent / "results" / "module03_structures.yaml"
    out.parent.mkdir(exist_ok=True)
    out.write_text(yaml.safe_dump(report, allow_unicode=True, sort_keys=False), encoding="utf-8")
    print(yaml.safe_dump(report["summary"], allow_unicode=True))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
