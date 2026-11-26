"""Analyse d'une séquence peptidique avec ses extrémités.

Notations acceptées :

- une lettre : `RGDS`, `Ac-RGDS-NH2`, `H-RGDS-OH` ;
- trois lettres séparées par des tirets : `Arg-Gly-Asp-Ser`, `Ac-Arg-Gly-Asp-Ser-NH2` ;
- notation complète : `NH2-RGDS-CO-NH2`, `H2N-Arg-Gly-CO-OH`.

Règle de lecture : un unique bloc central est lu en code à une lettre, sauf s'il
est écrit en casse « titre » (`Ala`), auquel cas il est lu comme un code à trois
lettres. `ALA` est donc lu A-L-A. Plusieurs blocs séparés par des tirets sont
toujours lus en codes à trois lettres. Sans extrémité déclarée, les valeurs par
défaut de la bibliothèque s'appliquent (amine libre, acide libre) et sont
signalées dans `Peptide.metadata`.
"""

from __future__ import annotations

from pepforge.core.peptide import Peptide, Residue
from pepforge.exceptions import SequenceParseError
from pepforge.residues.library import ResidueLibrary, load_default_library

_MAX_C_TOKENS = 3  # "CO-NH-CH3" est l'alias C-terminal le plus long


def _split_c_terminus(tokens: list[str], library: ResidueLibrary) -> tuple[list[str], str | None]:
    for size in range(min(_MAX_C_TOKENS, len(tokens) - 1), 0, -1):
        candidate = "-".join(tokens[-size:])
        if library.is_c_alias(candidate):
            return tokens[:-size], candidate
    return tokens, None


def _parse_body(tokens: list[str], library: ResidueLibrary) -> list[Residue]:
    if len(tokens) == 1:
        token = tokens[0]
        if len(token) == 3 and token.istitle() and library.has_code(token):
            return [library.by_code(token)]
        return [library.by_one_letter(letter) for letter in token]
    return [library.by_code(token) for token in tokens]


def parse_sequence(text: str, library: ResidueLibrary | None = None) -> Peptide:
    """Construit un `Peptide` à partir d'une séquence texte."""
    lib = library or load_default_library()
    cleaned = "".join(text.split())
    if not cleaned:
        raise SequenceParseError("Séquence vide")
    tokens = [t for t in cleaned.split("-")]
    if any(not t for t in tokens):
        raise SequenceParseError(f"Tiret mal placé dans la séquence : {text!r}")

    metadata: dict[str, str] = {"input": text}

    n_alias: str | None = None
    if len(tokens) > 1 and lib.is_n_alias(tokens[0]):
        n_alias, tokens = tokens[0], tokens[1:]

    c_alias: str | None = None
    if len(tokens) > 1:
        tokens, c_alias = _split_c_terminus(tokens, lib)

    if not tokens:
        raise SequenceParseError(f"Aucun résidu dans la séquence : {text!r}")

    if n_alias is None:
        metadata["n_terminus_source"] = "défaut de la bibliothèque"
    if c_alias is None:
        metadata["c_terminus_source"] = "défaut de la bibliothèque"

    return Peptide(
        residues=_parse_body(tokens, lib),
        n_terminus=lib.n_terminus(n_alias),
        c_terminus=lib.c_terminus(c_alias),
        metadata=metadata,
    )
