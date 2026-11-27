"""Analyse d'une séquence peptidique avec ses extrémités.

Notations acceptées :

- une lettre : `RGDS`, `Ac-RGDS-NH2`, `H-RGDS-OH` ;
- trois lettres séparées par des tirets : `Arg-Gly-Asp-Ser`, `Ac-Arg-Gly-Asp-Ser-NH2` ;
- notation complète : `NH2-RGDS-CO-NH2`, `H2N-Arg-Gly-CO-OH` ;
- résidu désigné par son code de bibliothèque, entre crochets, dans les deux
  notations : `Ac-RG[Aib]S-NH2`, `Ac-Arg-Gly-[Aib]-Ser-NH2`.

Règle de lecture, appliquée dans cet ordre :

1. Les espaces sont ignorés, les crochets sont vérifiés, puis la séquence est
   découpée en blocs aux tirets.
2. Le premier bloc est une extrémité N s'il en est un alias connu (`Ac`, `H`,
   `NH2`...). Les derniers blocs, jusqu'à trois pour `CO-NH-CH3`, forment une
   extrémité C s'ils en sont un alias connu.
3. Un unique bloc central est lu en codes à une lettre, sauf s'il est écrit en
   casse « titre » (`Ala`) et qu'il est un code de la bibliothèque, auquel cas il
   est lu comme un code à trois lettres. `ALA` est donc lu A-L-A. Comme `ALA` est
   aussi le code de l'alanine, ce cas est **signalé** : la lecture retenue et la
   façon d'écrire l'autre sont inscrites dans `Peptide.metadata`
   (`AMBIGUOUS_READING`), et le contrôle `peptide.unambiguous_reading` produit un
   avertissement. Si la lecture en une lettre échoue (`GLU`, où `U` n'est pas un
   code standard), l'erreur dit comment écrire le résidu seul (`Glu`).
4. Plusieurs blocs centraux sont toujours lus en codes à trois lettres.
5. Un code entre crochets est un code de la bibliothèque, standard ou non
   (`[Aib]`, `[Asp]`), sans distinction de casse. Dans un bloc à une lettre, il
   prend la place d'une lettre ; en notation à trois lettres, il occupe un bloc
   entier, et les crochets y sont facultatifs (`Arg-Aib-Ser` vaut
   `Arg-[Aib]-Ser`). Les crochets ne s'imbriquent pas, ne sont jamais vides et
   ne contiennent pas de tiret ; un bloc à trois lettres ne mêle pas lettres et
   crochets (`Arg[Aib]` est refusé). Ces crochets reprennent la convention de
   la notation HELM pour les monomères désignés par plusieurs caractères.
6. Sans extrémité déclarée, les valeurs par défaut de la bibliothèque
   s'appliquent (amine libre, acide libre) et sont signalées dans
   `Peptide.metadata`.

Un code entre crochets n'est reconnu que s'il figure dans la bibliothèque
chargée : la bibliothèque par défaut ne contient encore que les vingt acides
aminés standards, et un code absent lève `UnknownResidueError`.
"""

from __future__ import annotations

from pepforge.core.peptide import (
    AMBIGUOUS_READING,
    C_TERMINUS_DEFAULTED,
    N_TERMINUS_DEFAULTED,
    Peptide,
    Residue,
)
from pepforge.exceptions import SequenceParseError, UnknownResidueError
from pepforge.residues.library import ResidueLibrary, load_default_library

_MAX_C_TOKENS = 3  # "CO-NH-CH3" est l'alias C-terminal le plus long
_DEFAULT_SOURCE = "défaut de la bibliothèque"


def _check_brackets(cleaned: str) -> None:
    """Refuse les crochets mal formés avant tout découpage, avec un message précis.

    Les positions données sont celles des caractères de la séquence une fois
    les espaces retirés.
    """
    opened: int | None = None
    for index, char in enumerate(cleaned):
        if char == "[":
            if opened is not None:
                raise SequenceParseError(
                    f"Crochets imbriqués (caractère {index + 1}) dans {cleaned!r}"
                )
            opened = index
        elif char == "]":
            if opened is None:
                raise SequenceParseError(
                    f"Crochet fermant sans crochet ouvrant (caractère {index + 1}) dans {cleaned!r}"
                )
            if index == opened + 1:
                raise SequenceParseError(
                    f"Crochets vides (caractère {opened + 1}) dans {cleaned!r}"
                )
            opened = None
        elif char == "-" and opened is not None:
            raise SequenceParseError(
                f"Crochet ouvert (caractère {opened + 1}) et non refermé avant le tiret "
                f"(caractère {index + 1}) dans {cleaned!r} : un code entre crochets ne "
                "contient jamais de tiret"
            )
    if opened is not None:
        raise SequenceParseError(
            f"Crochet ouvert (caractère {opened + 1}) jamais fermé dans {cleaned!r}"
        )


def _split_c_terminus(tokens: list[str], library: ResidueLibrary) -> tuple[list[str], str | None]:
    for size in range(min(_MAX_C_TOKENS, len(tokens) - 1), 0, -1):
        candidate = "-".join(tokens[-size:])
        if library.is_c_alias(candidate):
            return tokens[:-size], candidate
    return tokens, None


def _one_letter_block(token: str, library: ResidueLibrary) -> list[Residue]:
    residues: list[Residue] = []
    index = 0
    while index < len(token):
        if token[index] == "[":
            end = token.index("]", index)  # présent : crochets déjà vérifiés
            residues.append(library.by_code(token[index + 1 : end]))
            index = end + 1
        else:
            residues.append(library.by_one_letter(token[index]))
            index += 1
    return residues


def _three_letter_token(token: str, library: ResidueLibrary) -> Residue:
    if "[" not in token and "]" not in token:
        return library.by_code(token)
    if token.startswith("[") and token.endswith("]") and token.count("[") == 1:
        return library.by_code(token[1:-1])
    raise SequenceParseError(
        f"Bloc {token!r} : en notation à trois lettres, un code entre crochets occupe "
        "un bloc entier (ex. Arg-[Aib]-Ser)"
    )


def _parse_body(tokens: list[str], library: ResidueLibrary) -> list[Residue]:
    if len(tokens) == 1:
        token = tokens[0]
        if "[" not in token and len(token) == 3 and token.istitle() and library.has_code(token):
            return [library.by_code(token)]
        return _one_letter_block(token, library)
    return [_three_letter_token(token, library) for token in tokens]


def _ambiguous_block(tokens: list[str], library: ResidueLibrary) -> str | None:
    """Bloc unique hors casse « titre » qui est aussi un code à trois lettres connu."""
    if len(tokens) != 1:
        return None
    token = tokens[0]
    if "[" in token or len(token) != 3 or token.istitle() or not library.has_code(token):
        return None
    return token


def parse_sequence(text: str, library: ResidueLibrary | None = None) -> Peptide:
    """Construit un `Peptide` à partir d'une séquence texte (règle de lecture ci-dessus)."""
    lib = library or load_default_library()
    cleaned = "".join(text.split())
    if not cleaned:
        raise SequenceParseError("Séquence vide")
    _check_brackets(cleaned)
    tokens = cleaned.split("-")
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
        metadata[N_TERMINUS_DEFAULTED] = _DEFAULT_SOURCE
    if c_alias is None:
        metadata[C_TERMINUS_DEFAULTED] = _DEFAULT_SOURCE

    ambiguous = _ambiguous_block(tokens, lib)
    try:
        residues = tuple(_parse_body(tokens, lib))
    except UnknownResidueError as error:
        if ambiguous is None:
            raise
        raise UnknownResidueError(
            f"{error}. « {ambiguous} » est lu en codes à une lettre ; pour le résidu "
            f"{ambiguous.upper()} seul, écrire « {ambiguous.title()} »"
        ) from None
    if ambiguous is not None:
        metadata[AMBIGUOUS_READING] = (
            f"« {ambiguous} » lu en codes à une lettre ({'-'.join(r.code for r in residues)}) ; "
            f"pour le résidu {ambiguous.upper()} seul, écrire « {ambiguous.title()} »"
        )

    return Peptide(
        residues=residues,
        n_terminus=lib.n_terminus(n_alias),
        c_terminus=lib.c_terminus(c_alias),
        metadata=metadata,
    )
