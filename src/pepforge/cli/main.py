"""Point d'entrée de la commande `pepforge`."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

import yaml

from pepforge import __version__
from pepforge.config import PrepConfig
from pepforge.core import C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED, Peptide
from pepforge.exceptions import (
    ConfigurationError,
    ExtractionRefusedError,
    PepForgeError,
    SequenceParseError,
)
from pepforge.io.peptide_extraction import (
    BREAK_DECISIONS,
    COMPONENT_DECISIONS,
    CONTACT_DECISIONS,
    STRUCTURE_SEGMENT,
    UNOBSERVED_TERMINUS_DECISIONS,
    ExtractionDecisions,
    ExtractionReport,
    extract_peptides,
)
from pepforge.io.retrieval import DownloadLimits
from pepforge.io.sequence import SequenceNotation, format_sequence, parse_sequence
from pepforge.io.structure import read_structure, read_structure_url
from pepforge.io.structure_inventory import inventory
from pepforge.report.journal import Journal

_NOTATION_CHOICES: dict[str, SequenceNotation] = {
    "1": SequenceNotation.ONE_LETTER,
    "3": SequenceNotation.THREE_LETTER,
    SequenceNotation.ONE_LETTER.value: SequenceNotation.ONE_LETTER,
    SequenceNotation.THREE_LETTER.value: SequenceNotation.THREE_LETTER,
}
_NOTATION_HELP = (
    "mode de lecture, à préciser dès le départ : 1 = code à une lettre (monolitère), "
    "3 = code à trois lettres (trilitère)"
)


def _notation(value: str | None) -> SequenceNotation:
    """Mode de lecture demandé en ligne de commande ; jamais deviné."""
    if value is None:
        raise SequenceParseError(
            "Mode de lecture non précisé : ajouter --notation 1 (code à une lettre, "
            "monolitère) ou --notation 3 (code à trois lettres, trilitère)"
        )
    try:
        return _NOTATION_CHOICES[value]
    except KeyError:
        raise SequenceParseError(
            f"Mode de lecture inconnu : {value!r} ; --notation 1 (code à une lettre, "
            "monolitère) ou --notation 3 (code à trois lettres, trilitère)"
        ) from None


def _termini_warnings(peptide: Peptide) -> list[str]:
    return [
        f"Attention : extrémité {side} non déclarée, valeur par défaut."
        for side, key in (("N", N_TERMINUS_DEFAULTED), ("C", C_TERMINUS_DEFAULTED))
        if key in peptide.metadata
    ]


def _cmd_parse(args: argparse.Namespace) -> int:
    notation = _notation(args.notation)
    peptide = parse_sequence(args.sequence, notation=notation)
    one_letter = format_sequence(peptide, SequenceNotation.ONE_LETTER)
    three_letter = format_sequence(peptide, SequenceNotation.THREE_LETTER)
    print(f"Mode de lecture      : {notation.label}")
    print(f"Code à une lettre    : {one_letter}")
    print(f"Code à trois lettres : {three_letter}")
    print(f"Longueur             : {len(peptide)}")
    print(f"Extrémité N          : {peptide.n_terminus.label}")
    print(f"Extrémité C          : {peptide.c_terminus.label}")
    composition = ", ".join(f"{k}×{v}" for k, v in sorted(peptide.composition().items()))
    print(f"Composition          : {composition}")
    for warning in _termini_warnings(peptide):
        print(warning)
    return 0


def _cmd_convert(args: argparse.Namespace) -> int:
    """Bascule la séquence dans l'autre code ; seul le résultat va sur la sortie standard."""
    notation = _notation(args.notation)
    peptide = parse_sequence(args.sequence, notation=notation)
    print(format_sequence(peptide, notation.other))
    for warning in _termini_warnings(peptide):
        print(warning, file=sys.stderr)
    return 0


def _extraction_summary(report: ExtractionReport) -> dict[str, object]:
    """Rapport d'extraction lisible : une ligne par jonction et par composant (eaux comptées)."""
    data = report.to_dict()
    junctions = []
    for junction in report.junctions:
        atoms = "; ".join(
            f"{m.c_atom}–{m.n_atom} {m.distance:.3f} Å" for m in junction.measurements
        )
        junctions.append(
            f"{junction.before} → {junction.after} | {atoms or 'non mesurable'} | "
            f"{junction.category.value} | {junction.decision.value}"
        )
    waters = [c for c in report.components if c.category.value == "water"]
    return {
        "chaîne": report.chain_id,
        "acceptée": report.accepted,
        "critère": data["criterion"],
        "décisions": data["decisions"],
        "jonctions_par_catégorie": data["junctions_by_category"],
        "jonctions": junctions,
        "composants": [
            f"{c.residue} | {c.category.value} | {c.decision.value} | {c.evidence}"
            + (
                f" | distance minimale à la chaîne {c.min_distance_to_chain:.3f} Å"
                if c.min_distance_to_chain is not None
                else ""
            )
            + (f" | {'; '.join(c.bonds)}" if c.bonds else "")
            for c in report.components
            if c not in waters
        ],
        "eaux": len(waters),
        "liaisons": [
            f"{b.first.label()} – {b.second.label()} | "
            + (f"{b.distance:.3f} Å" if b.distance is not None else "non mesurée")
            + f" | {b.evidence.value}"
            + (f" ({'+'.join(b.records)})" if b.records else "")
            + f" | {b.certainty.value} | {b.status.value} | {b.reason}"
            for b in report.bonds
        ],
        "critère_covalent": data["covalent_criteria"],
        "extrémités_non_observées": data["unobserved_termini"],
        "résidus_déclarés_non_observés": len(report.declared_missing_residues),
        "atomes_déclarés_non_observés": list(report.declared_missing_atoms),
        "atomes_de_squelette_absents": list(report.absent_backbone_atoms),
        "segments": data["segments"],
        "refus": list(report.refusals),
    }


def _named(values: Sequence[str] | None, option: str) -> dict[str, str]:
    """`--component A:DOT201=exclude_as_unbound` répété → tableau ; un doublon est refusé."""
    named: dict[str, str] = {}
    for value in values or ():
        label, sep, decision = value.rpartition("=")
        if not sep or not label or not decision:
            raise ConfigurationError(f"{option} {value!r} : forme attendue ÉTIQUETTE=DÉCISION")
        if label in named:
            raise ConfigurationError(f"{option} : {label} nommé deux fois")
        named[label] = decision
    return named


def _cmd_inspect(args: argparse.Namespace) -> int:
    """Importe une structure (fichier ou URL), sans la modifier, et affiche son inventaire."""
    journal = Journal()
    if args.url:
        limits = DownloadLimits(
            timeout_s=args.timeout,
            max_bytes=int(args.max_mb * 1024 * 1024),
            allow_http=args.allow_http,
        )
        structure = read_structure_url(
            args.url, journal=journal, model=args.model, fmt=args.format, limits=limits
        )
    else:
        if not args.path:
            raise SequenceParseError("Indiquer un fichier, ou une adresse avec --url")
        structure = read_structure(args.path, journal=journal, model=args.model, fmt=args.format)
    report: dict[str, object] = {
        "provenance": dict(journal.entries[0].parameters),
        "inventaire": inventory(structure).to_dict(),
    }
    if args.chain:
        decisions = ExtractionDecisions.from_dict(
            {
                "breaks": args.breaks,
                "unobserved_termini": args.unobserved_termini,
                "components": _named(args.component, "--component"),
                "contacts": _named(args.contact, "--contact"),
            },
            "décisions de la ligne de commande",
        )
        try:
            result = extract_peptides(structure, args.chain, journal=journal, decisions=decisions)
        except ExtractionRefusedError as exc:
            if isinstance(exc.report, ExtractionReport):
                report["extraction"] = _extraction_summary(exc.report)
                print(
                    yaml.safe_dump(report, allow_unicode=True, sort_keys=False, width=1000), end=""
                )
            raise
        report["extraction"] = _extraction_summary(result.report)
        report["peptides"] = [
            {
                "segment": peptide.metadata.get(STRUCTURE_SEGMENT),
                "résidus": "-".join(residue.code for residue in peptide.residues),
                "extrémité_N": peptide.n_terminus.label,
                "extrémité_C": peptide.c_terminus.label,
                "ponts_disulfure": [list(link.positions) for link in peptide.links],
                "avertissements": _termini_warnings(peptide),
            }
            for peptide in result.peptides
        ]
    print(yaml.safe_dump(report, allow_unicode=True, sort_keys=False, width=1000), end="")
    return 0


def _cmd_config_template(_: argparse.Namespace) -> int:
    config = PrepConfig()
    print("# Modèle de configuration PepForge : toutes les décisions sont à prendre.")
    print(yaml.safe_dump(config.to_dict(), allow_unicode=True, sort_keys=False), end="")
    print("# À décider : " + "; ".join(config.missing_decisions()))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pepforge",
        description="PepForge — Universal Peptide Design, Engineering and Simulation Platform",
    )
    parser.add_argument("--version", action="version", version=f"pepforge {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_parse = sub.add_parser("parse", help="Analyser une séquence et ses extrémités")
    p_parse.add_argument(
        "sequence", help="ex. Ac-RGDS-NH2 (code à une lettre) ou Ac-Arg-Gly-Asp-Ser-NH2"
    )
    p_parse.add_argument("-n", "--notation", help=_NOTATION_HELP)
    p_parse.set_defaults(func=_cmd_parse)

    p_convert = sub.add_parser(
        "convert", help="Réécrire une séquence dans l'autre code (une lettre ↔ trois lettres)"
    )
    p_convert.add_argument("sequence", help="séquence à réécrire, dans le code donné par -n")
    p_convert.add_argument("-n", "--notation", help=_NOTATION_HELP)
    p_convert.set_defaults(func=_cmd_convert)

    p_inspect = sub.add_parser(
        "inspect",
        help="Importer une structure (fichier local ou URL) et afficher son inventaire",
    )
    p_inspect.add_argument("path", nargs="?", help="fichier PDB ou mmCIF (éventuellement .gz)")
    p_inspect.add_argument("--url", help="adresse du fichier, téléchargé sous contrôle")
    p_inspect.add_argument(
        "--model", type=int, help="modèle à lire (obligatoire si le fichier en a plusieurs)"
    )
    p_inspect.add_argument(
        "--format", choices=["pdb", "mmcif"], help="format attendu, vérifié sur le contenu"
    )
    p_inspect.add_argument("--chain", help="établir le peptide décrit par cette chaîne")
    p_inspect.add_argument(
        "--breaks",
        choices=[d.value for d in BREAK_DECISIONS],
        default="reject",
        help="ruptures de chaîne : refuser (défaut) ou découper en segments",
    )
    p_inspect.add_argument(
        "--unobserved-termini",
        choices=[d.value for d in UNOBSERVED_TERMINUS_DECISIONS],
        default="reject",
        help="résidus non observés aux extrémités : refuser (défaut) ou garder le segment observé",
    )
    p_inspect.add_argument(
        "--component",
        action="append",
        metavar="ÉTIQUETTE=DÉCISION",
        help="décision pour un composant nommé, par exemple A:DOT201=exclude_as_unbound "
        f"(décisions : {', '.join(d.value for d in COMPONENT_DECISIONS)}) ; répétable",
    )
    p_inspect.add_argument(
        "--contact",
        action="append",
        metavar="PAIRE=DÉCISION",
        help="décision pour un contact covalent possible, par exemple "
        f"'A:CYS3–A:CYS10=disulfide' (décisions : {', '.join(d.value for d in CONTACT_DECISIONS)}) "
        "; répétable",
    )
    p_inspect.add_argument(
        "--timeout", type=float, default=DownloadLimits.timeout_s, help="délai réseau (s)"
    )
    p_inspect.add_argument(
        "--max-mb",
        type=float,
        default=DownloadLimits.max_bytes / (1024 * 1024),
        help="taille maximale acceptée (Mio)",
    )
    p_inspect.add_argument(
        "--allow-http", action="store_true", help="accepter une adresse http non chiffrée"
    )
    p_inspect.set_defaults(func=_cmd_inspect)

    p_tpl = sub.add_parser("config-template", help="Afficher un modèle de configuration")
    p_tpl.set_defaults(func=_cmd_config_template)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        code: int = args.func(args)
    except PepForgeError as exc:
        print(f"Erreur : {exc}", file=sys.stderr)
        return 1
    return code


if __name__ == "__main__":
    sys.exit(main())
