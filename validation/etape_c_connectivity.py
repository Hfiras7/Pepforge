"""Validation de l'étape C (liaisons covalentes) sur les structures réelles du module 3.

Ce script n'est pas un test pytest : il demande le réseau. Il se relance tel quel :

    pip install -e ".[validation]"
    python validation/etape_c_connectivity.py

Mêmes fichiers, mêmes adresses figées que la validation du module 3
(`validation/module03_structures.py` : suite de tests de Biopython 1.88, copies
d'entrées de la PDB). Les résultats du module 3 servent de référence :
`validation/results/module03_structures.yaml`, produit au commit `1cd2cab`.

Contrôles, pour chaque chaîne peptidique déjà examinée au module 3 :

1. **Résultat du module 3 conservé** : chaque motif de refus du module 3 est
   toujours là, et aucun composant ne change de catégorie.
2. **Nouveaux refus expliqués** : chaque motif nouveau figure dans
   `EXPECTED_NEW`, table établie à la main après examen de chaque cas ; un
   motif nouveau absent de la table fait échouer la validation.
3. **Pas de faux positif systématique** : toutes les liaisons possibles (sans
   déclaration) trouvées par la géométrie, sur tous les fichiers, sont listées ;
   chacune doit être dans `EXPECTED_CONTACTS`.
4. **Liaisons déclarées confirmées** : chaque liaison covalente ou disulfure
   déclarée qui touche une chaîne peptidique est confrontée aux coordonnées.
5. **Ligands non liés et ions** : ceux du module 3 le restent.

Puis des **contrôles positifs** sur fichiers réels modifiés (on retire une
information du fichier et on vérifie que le composant n'est pas perdu) :

6. 7DDO sans LINK (CONECT gardés) : les glycanes restent rattachés, par CONECT ;
7. 7DDO sans LINK ni CONECT : les glycanes deviennent `possibly_bound`, refusés ;
   la décision `exclude_as_unbound` les laisse hors du peptide, et le dit ;
8. 6WQA avec la décision `disulfide` : les quatre ponts entrent dans le peptide ;
9. 4ZHL chaîne U avec la décision `disulfide` : l'extraction redevient acceptée,
   la séquence lue est celle de `_entity_poly` (contrôle que le module 3 faisait
   avant que l'étape C refuse cette chaîne sans décision), et les deux ponts
   s'ajoutent aux trois ponts de la chaîne U déclarés par le fichier (le
   quatrième pont déclaré, CYS1–CYS10, appartient à la chaîne P).

Les résultats sont écrits dans `validation/results/etape_c_connectivity.yaml`.
"""

from __future__ import annotations

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import yaml

sys.path.insert(0, str(Path(__file__).parent))

from module03_structures import BASE, FILES, _entity_sequence

from pepforge.exceptions import ExtractionRefusedError, PepForgeError, StructureError
from pepforge.io import (
    BondEvidence,
    BondStatus,
    ExtractionDecisions,
    ExtractionReport,
    analyze_chain,
    extract_peptides,
    read_structure,
)
from pepforge.io.peptide_extraction import STRUCTURE_SET_ASIDE
from pepforge.io.retrieval import download_file
from pepforge.report.journal import Journal

RESULTS = Path(__file__).parent / "results"
MODULE03 = RESULTS / "module03_structures.yaml"

# Motifs de refus nouveaux, établis après examen de chaque cas (bilan de l'étape C) :
EXPECTED_NEW: dict[tuple[str, str], list[str]] = {
    # Quatre contacts SG–SG à 2,02–2,03 Å, sans struct_conn dans cette copie du fichier :
    # des ponts disulfure que le module 3 ignorait en silence.
    ("6WQA.cif", "A"): ["contact_in_chain"] * 4,
    # Deux contacts SG–SG à 1,51 et 1,64 Å, sans struct_conn, alors que les trois autres
    # ponts de la chaîne sont déclarés : à trancher (pont non déclaré ou géométrie fautive).
    ("4ZHL.cif", "U"): ["contact_in_chain"] * 2,
    # Huit CONECT de cette copie désignent des numéros de série absents (les atomes de
    # MSE151 ont été renumérotés dans la copie de Biopython) : connectivité incertaine.
    ("1A8O.pdb", "A"): ["conect_unresolved"],
}

# Toutes les liaisons possibles attendues (géométrie seule), avec la distance arrondie :
EXPECTED_CONTACTS: dict[str, set[str]] = {
    "6WQA.cif": {
        "A:CYS71:SG–A:CYS159:SG",
        "A:CYS74:SG–A:CYS146:SG",
        "A:CYS77:SG–A:CYS166:SG",
        "A:CYS259:SG–A:CYS262:SG",
    },
    "4ZHL.cif": {"U:CYS50:SG–U:CYS111:SG", "U:CYS136:SG–U:CYS201:SG"},
}

WQA_DISULFIDES = {
    "A:CYS71–A:CYS159": "disulfide",
    "A:CYS74–A:CYS146": "disulfide",
    "A:CYS77–A:CYS166": "disulfide",
    "A:CYS259–A:CYS262": "disulfide",
}

ZHL_DISULFIDES = {"U:CYS50–U:CYS111": "disulfide", "U:CYS136–U:CYS201": "disulfide"}


def _new_causes(report: ExtractionReport) -> list[str]:
    """Motifs de refus propres à l'étape C, par nature."""
    causes: list[str] = []
    for reason in report.refusals:
        if reason.startswith("contact covalent possible dans la chaîne"):
            causes.append("contact_in_chain")
        elif reason.startswith("contact covalent possible avec une autre chaîne"):
            causes.append("contact_between_chains")
        elif reason.startswith("connectivité du fichier incertaine"):
            causes.append("conect_unresolved")
        elif reason.startswith("composant ") and "au contact covalent" in reason:
            causes.append("possibly_bound")
        elif "contredite par les coordonnées" in reason:
            causes.append("declared_contradicted")
    return causes


def _check(results: dict[str, Any], name: str, ok: bool, detail: str) -> None:
    results["checks"].append({"check": name, "ok": ok, "detail": detail})


def _module03() -> dict[str, dict[str, Any]]:
    data = yaml.safe_load(MODULE03.read_text(encoding="utf-8"))
    return {entry["file"]: entry for entry in data["files"]}


def _download(name: str, workdir: Path) -> Path:
    retrieved = download_file(BASE + name)
    path = workdir / name
    path.write_bytes(retrieved.content)
    return path


def _read(path: Path) -> Any:
    try:
        return read_structure(path, journal=Journal())
    except StructureError as exc:
        if "modèles" in str(exc):
            return read_structure(path, journal=Journal(), model=1)
        raise


def validate_file(name: str, path: Path, reference: dict[str, Any]) -> dict[str, Any]:
    results: dict[str, Any] = {"file": name, "checks": [], "chains": {}}
    try:
        structure = _read(path)
    except PepForgeError as exc:
        expected_unreadable = not reference.get("peptides")
        _check(
            results,
            "lecture",
            expected_unreadable,
            f"lecture refusée ({str(exc)[:160]}) ; module 3 : "
            + ("refusée aussi" if expected_unreadable else "lue"),
        )
        return results
    contacts_found: set[str] = set()
    for chain, old in (reference.get("peptides") or {}).items():
        if isinstance(old, str):
            # Chaîne refusée d'emblée au module 3 (non peptidique) : doit l'être encore.
            try:
                analyze_chain(structure, chain)
                again = "acceptée à l'analyse"
            except StructureError as exc:
                again = f"refusé ({type(exc).__name__}) : {exc}"
            _check(results, f"{chain} : refus du module 3 conservé", again == old, again[:120])
            continue
        report = analyze_chain(structure, chain)
        new = _new_causes(report)
        old_causes = old.get("refusal_causes", [])
        summary: dict[str, Any] = {
            "accepted": report.accepted,
            "new_refusal_causes": new,
            "components": {
                c.residue: c.category.value
                for c in report.components
                if c.category.value != "water"
            },
            "bonds": [
                f"{b.label()} | "
                + (f"{b.distance:.3f} Å" if b.distance is not None else "non mesurée")
                + f" | {b.evidence.value} {'+'.join(b.records)} | {b.certainty.value}"
                for b in report.bonds
            ],
        }
        results["chains"][chain] = summary
        # 1. Module 3 conservé.
        changed = {
            residue: (category, summary["components"].get(residue))
            for residue, category in (old.get("components") or {}).items()
            if summary["components"].get(residue) != category
        }
        module3_refused = bool(old_causes)
        still = (not module3_refused) or (not report.accepted)
        _check(
            results,
            f"{chain} : résultat du module 3 conservé",
            not changed and still,
            f"motifs du module 3 {old_causes} ; catégories changées : {changed or 'aucune'}",
        )
        # 2. Nouveaux refus expliqués.
        expected = EXPECTED_NEW.get((name, chain), [])
        _check(
            results,
            f"{chain} : nouveaux refus expliqués",
            sorted(new) == sorted(expected),
            f"trouvés {new} ; attendus {expected}",
        )
        # 3. Liaisons possibles.
        for bond in report.bonds:
            if bond.evidence is BondEvidence.GEOMETRIC:
                contacts_found.add(bond.label())
        # 4. Liaisons déclarées confirmées.
        declared = [b for b in report.bonds if b.evidence is BondEvidence.DECLARED]
        covalent = [b for b in declared if b.kind in ("covalent", "disulfide", "unknown")]
        bad = [b.label() for b in covalent if b.status is not BondStatus.ACCEPTED]
        _check(
            results,
            f"{chain} : liaisons covalentes déclarées confirmées par les coordonnées",
            not bad,
            f"{len(covalent) - len(bad)}/{len(covalent)} confirmées"
            + (f" ; non confirmées : {bad}" if bad else ""),
        )
        # 5. Ligands non liés et ions.
        unbound = {
            r: c for r, c in (old.get("components") or {}).items() if c in ("unbound_ligand", "ion")
        }
        kept = all(summary["components"].get(r) == c for r, c in unbound.items())
        _check(
            results,
            f"{chain} : ligands non liés et ions inchangés",
            kept,
            f"{len(unbound)} composant(s) : {sorted(unbound) or 'aucun'}",
        )
    expected_contacts = EXPECTED_CONTACTS.get(name, set())
    _check(
        results,
        "liaisons possibles (géométrie seule) : toutes attendues",
        contacts_found == expected_contacts,
        f"trouvées {sorted(contacts_found)} ; attendues {sorted(expected_contacts)}",
    )
    return results


def _without(path: Path, records: tuple[str, ...], target: Path) -> Path:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    target.write_text(
        "".join(line for line in lines if not line.startswith(records)), encoding="utf-8"
    )
    return target


def positive_controls(workdir: Path) -> dict[str, Any]:
    results: dict[str, Any] = {"file": "contrôles positifs", "checks": []}
    ddo = workdir / "7DDO.pdb"
    glycans = {"A:NAG902", "A:NAG903", "A:NAG904"}

    # 6. Sans LINK, CONECT gardés.
    no_link = read_structure(
        _without(ddo, ("LINK  ",), workdir / "7DDO_sans_LINK.pdb"), journal=Journal()
    )
    report = analyze_chain(no_link, "A")
    found = {c.residue: c.category.value for c in report.components if c.residue in glycans}
    records = {
        b.label(): "+".join(b.records)
        for b in report.bonds
        if glycans & {b.first.residue.label(), b.second.residue.label()}
    }
    _check(
        results,
        "7DDO sans LINK : glycanes rattachés par CONECT",
        set(found.values()) == {"side_chain_modification"}
        and len(found) == 3
        and set(records.values()) == {"CONECT"}
        and not report.accepted,
        f"catégories {found} ; sources {records}",
    )

    # 7. Sans LINK ni CONECT.
    bare = read_structure(
        _without(ddo, ("LINK  ", "CONECT"), workdir / "7DDO_sans_LINK_CONECT.pdb"),
        journal=Journal(),
    )
    report = analyze_chain(bare, "A")
    found = {c.residue: c.category.value for c in report.components if c.residue in glycans}
    contacts = [
        f"{b.label()} {b.distance:.2f} Å"
        for b in report.bonds
        if b.evidence is BondEvidence.GEOMETRIC and b.distance is not None
    ]
    _check(
        results,
        "7DDO sans LINK ni CONECT : glycanes possibly_bound, extraction refusée",
        set(found.values()) == {"possibly_bound"} and len(found) == 3 and not report.accepted,
        f"catégories {found} ; contacts {contacts}",
    )
    zinc = [c for c in report.components if c.name == "ZN"]
    _check(
        results,
        "7DDO sans LINK ni CONECT : le zinc reste un ion",
        [c.category.value for c in zinc] == ["ion"],
        f"{[(c.residue, c.category.value) for c in zinc]}",
    )
    decisions = ExtractionDecisions(components={g: "exclude_as_unbound" for g in glycans})
    report = analyze_chain(bare, "A", decisions=decisions)
    decided = {c.residue: c.decision.value for c in report.components if c.residue in glycans}
    _check(
        results,
        "7DDO sans LINK ni CONECT : exclude_as_unbound appliqué et rapporté",
        set(decided.values()) == {"exclude_as_unbound"} and report.accepted,
        f"décisions {decided} ; refus restants {list(report.refusals)}",
    )
    result = extract_peptides(bare, "A", journal=Journal(), decisions=decisions)
    marked = result.peptides[0].metadata.get(STRUCTURE_SET_ASIDE, "")
    _check(
        results,
        "7DDO : composants écartés inscrits dans le peptide",
        all(g in marked for g in glycans),
        marked,
    )

    # 8. 6WQA avec les ponts admis.
    wqa = read_structure(workdir / "6WQA.cif", journal=Journal())
    try:
        extract_peptides(wqa, "A", journal=Journal())
        refused = False
    except ExtractionRefusedError:
        refused = True
    _check(results, "6WQA sans décision : refusée", refused, "contacts SG–SG non déclarés")
    split = extract_peptides(
        wqa,
        "A",
        journal=Journal(),
        decisions=ExtractionDecisions(breaks="split", contacts=WQA_DISULFIDES),
    )
    links = [link.describe() for p in split.peptides for link in p.links]
    _check(
        results,
        "6WQA avec la décision disulfide : quatre ponts dans les peptides",
        len(links) == 4,
        f"{links}",
    )

    # 9. 4ZHL chaîne U avec les deux ponts admis.
    zhl_path = workdir / "4ZHL.cif"
    zhl = read_structure(zhl_path, journal=Journal())
    chain_u = extract_peptides(
        zhl, "U", journal=Journal(), decisions=ExtractionDecisions(contacts=ZHL_DISULFIDES)
    )
    (peptide,) = chain_u.peptides
    read = "".join(r.one_letter or "?" for r in peptide.residues)
    declared = _entity_sequence(zhl_path, "U")
    _check(
        results,
        "4ZHL U avec la décision disulfide : séquence lue = séquence déclarée (_entity_poly)",
        chain_u.report.accepted and read == declared,
        f"{len(read)} résidus lus, {len(declared or '')} déclarés",
    )
    _check(
        results,
        "4ZHL U : trois ponts déclarés (module 3) + deux admis par décision",
        {link.positions for link in peptide.links}
        == {(31, 47), (167, 183), (194, 222), (39, 110), (135, 204)},
        f"{[link.describe() for link in peptide.links]}",
    )
    return results


def main() -> int:
    reference = _module03()
    out: list[dict[str, Any]] = []
    with TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        for name in FILES:
            path = _download(name, workdir)
            out.append(validate_file(name, path, reference.get(name, {})))
        out.append(positive_controls(workdir))
    checks = [c for entry in out for c in entry["checks"]]
    failed = [
        f"{entry['file']} : {c['check']}" for entry in out for c in entry["checks"] if not c["ok"]
    ]
    summary = {"checks": len(checks), "failed": len(failed), "failures": failed}
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "etape_c_connectivity.yaml").write_text(
        yaml.safe_dump(
            {"source": BASE, "reference": str(MODULE03.name), "files": out, "summary": summary},
            allow_unicode=True,
            sort_keys=False,
            width=1000,
        ),
        encoding="utf-8",
    )
    print(yaml.safe_dump(summary, allow_unicode=True, sort_keys=False), end="")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
