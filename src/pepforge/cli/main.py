"""Point d'entrée de la commande `pepforge`."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

import yaml

from pepforge import __version__
from pepforge.config import PrepConfig
from pepforge.core import C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED, Peptide
from pepforge.exceptions import PepForgeError, SequenceParseError
from pepforge.io.retrieval import DownloadLimits
from pepforge.io.sequence import SequenceNotation, format_sequence, parse_sequence
from pepforge.io.structure import read_structure, read_structure_url
from pepforge.io.structure_inventory import inventory, structure_to_peptide
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
        peptide = structure_to_peptide(structure, args.chain, journal=journal)
        report["peptide"] = {
            "chaîne": args.chain,
            "résidus": "-".join(residue.code for residue in peptide.residues),
            "extrémité_N": peptide.n_terminus.label,
            "extrémité_C": peptide.c_terminus.label,
            "ponts_disulfure": [list(link.positions) for link in peptide.links],
            "avertissements": _termini_warnings(peptide),
        }
    print(yaml.safe_dump(report, allow_unicode=True, sort_keys=False), end="")
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
