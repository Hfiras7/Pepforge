"""Point d'entrée de la commande `pepforge`."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

import yaml

from pepforge import __version__
from pepforge.config import PrepConfig
from pepforge.core import AMBIGUOUS_READING, C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED
from pepforge.exceptions import PepForgeError
from pepforge.io.sequence import parse_sequence


def _cmd_parse(args: argparse.Namespace) -> int:
    peptide = parse_sequence(args.sequence)
    print(f"Séquence (1 lettre)  : {peptide.one_letter_sequence}")
    print(f"Séquence (3 lettres) : {peptide.three_letter_sequence}")
    print(f"Longueur             : {len(peptide)}")
    print(f"Extrémité N          : {peptide.n_terminus.label}")
    print(f"Extrémité C          : {peptide.c_terminus.label}")
    composition = ", ".join(f"{k}×{v}" for k, v in sorted(peptide.composition().items()))
    print(f"Composition          : {composition}")
    for side, key in (("N", N_TERMINUS_DEFAULTED), ("C", C_TERMINUS_DEFAULTED)):
        if key in peptide.metadata:
            print(f"Attention : extrémité {side} non déclarée, valeur par défaut.")
    if AMBIGUOUS_READING in peptide.metadata:
        print(f"Attention : lecture à confirmer, {peptide.metadata[AMBIGUOUS_READING]}.")
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
    p_parse.add_argument("sequence", help="ex. Ac-RGDS-NH2 ou NH2-Arg-Gly-Asp-CO-NH2")
    p_parse.set_defaults(func=_cmd_parse)

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
