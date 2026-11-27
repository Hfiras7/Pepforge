"""Lecture et écriture d'une séquence peptidique avec ses extrémités.

Le **mode de lecture se déclare au départ** (revue du module 1) : code à une
lettre, dit monolitère (`SequenceNotation.ONE_LETTER`), ou code à trois
lettres, dit trilitère (`SequenceNotation.THREE_LETTER`). Rien n'est deviné :
`ALA` déclaré en code à une lettre est le tripeptide Ala-Leu-Ala, déclaré en
code à trois lettres c'est l'alanine seule.

Notations acceptées :

- code à une lettre : `RGDS`, `Ac-RGDS-NH2`, `H-RGDS-OH`, `NH2-RGDS-CO-NH2` ;
- code à trois lettres : `Arg-Gly-Asp-Ser`, `Ac-Arg-Gly-Asp-Ser-NH2`,
  `H2N-Arg-Gly-CO-OH` ;
- dans les deux codes, un résidu désigné par son code de bibliothèque entre
  crochets : `Ac-RG[Aib]S-NH2`, `Ac-Arg-Gly-[Aib]-Ser-NH2`.

Règle de lecture, appliquée dans cet ordre :

1. Les espaces sont ignorés, les crochets sont vérifiés, puis la séquence est
   découpée en blocs aux tirets.
2. Le premier bloc est une extrémité N s'il en est un alias connu (`Ac`, `H`,
   `NH2`...). Les derniers blocs, jusqu'à trois pour `CO-NH-CH3`, forment une
   extrémité C s'ils en sont un alias connu.
3. En code à une lettre, les résidus forment **un seul bloc**, sans tiret, et
   chaque lettre majuscule est un résidu. Une **minuscule est refusée** : dans
   l'usage courant, elle désigne un acide aminé D (le `f` de `c(RGDfK)` est la
   D-phénylalanine), que la bibliothèque ne décrit pas encore ; la lire comme
   un acide aminé L serait une erreur silencieuse.
4. En code à trois lettres, **chaque résidu occupe un bloc**, sans distinction
   de casse (`Arg`, `ARG`, `arg`).
5. Un code entre crochets est un code de la bibliothèque, standard ou non
   (`[Aib]`, `[Asp]`), sans distinction de casse. En code à une lettre, il
   prend la place d'une lettre ; en code à trois lettres, il occupe un bloc
   entier, et les crochets y sont facultatifs (`Arg-Aib-Ser` vaut
   `Arg-[Aib]-Ser`). Les crochets ne s'imbriquent pas, ne sont jamais vides et
   ne contiennent pas de tiret ; un bloc à trois lettres ne mêle pas lettres et
   crochets (`Arg[Aib]` est refusé). Ces crochets reprennent la convention de
   la notation HELM pour les monomères désignés par plusieurs caractères.
6. Sans extrémité déclarée, les valeurs par défaut de la bibliothèque
   s'appliquent (amine libre, acide libre) et sont signalées dans
   `Peptide.metadata`.

Le mode déclaré est inscrit dans `Peptide.metadata` (`SEQUENCE_NOTATION`).

**Conversion entre les deux codes** (`convert_sequence`), à la manière du
réglage `seq_view_format` de PyMOL, qui bascule l'affichage d'une séquence
entre codes à une et à trois lettres : la séquence est lue dans son code, puis
réécrite dans l'autre (`format_sequence`). L'écriture est canonique : résidus
standards selon l'IUPAC (`R` ou `Arg`), non standards entre crochets avec leur
code de bibliothèque (`[AIB]`), extrémités sous leur forme usuelle (le premier
alias de `termini.yaml`). Une extrémité non déclarée le reste. Toute écriture
est **relue et comparée** au peptide d'origine : une conversion qui changerait
la molécule est refusée, jamais rendue.

Un code entre crochets n'est reconnu que s'il figure dans la bibliothèque
chargée : la bibliothèque par défaut ne contient encore que les vingt acides
aminés standards, et un code absent lève `UnknownResidueError`.
"""

from __future__ import annotations

from enum import StrEnum

from pepforge.core.peptide import C_TERMINUS_DEFAULTED, N_TERMINUS_DEFAULTED, Peptide, Residue
from pepforge.exceptions import (
    PepForgeError,
    SequenceParseError,
    UnknownResidueError,
    UnsupportedFeatureError,
)
from pepforge.residues.library import ResidueLibrary, load_default_library

SEQUENCE_NOTATION = "sequence_notation"
"""Clé de `Peptide.metadata` : mode de lecture déclaré pour la séquence saisie."""

_MAX_C_TOKENS = 3  # "CO-NH-CH3" est l'alias C-terminal le plus long
_DEFAULT_SOURCE = "défaut de la bibliothèque"


class SequenceNotation(StrEnum):
    """Mode de lecture d'une séquence, déclaré au départ : rien n'est deviné."""

    ONE_LETTER = "one_letter"
    """Code à une lettre (monolitère) : `Ac-RGDS-NH2`."""

    THREE_LETTER = "three_letter"
    """Code à trois lettres (trilitère) : `Ac-Arg-Gly-Asp-Ser-NH2`."""

    @property
    def label(self) -> str:
        """Libellé lisible en français."""
        return _NOTATION_LABELS[self]

    @property
    def other(self) -> SequenceNotation:
        """L'autre code : celui vers lequel bascule une conversion."""
        if self is SequenceNotation.ONE_LETTER:
            return SequenceNotation.THREE_LETTER
        return SequenceNotation.ONE_LETTER


_NOTATION_LABELS: dict[SequenceNotation, str] = {
    SequenceNotation.ONE_LETTER: "code à une lettre (monolitère)",
    SequenceNotation.THREE_LETTER: "code à trois lettres (trilitère)",
}


def _as_notation(value: SequenceNotation | str) -> SequenceNotation:
    try:
        return SequenceNotation(value)
    except ValueError:
        accepted = " ; ".join(f"{n.value} = {n.label}" for n in SequenceNotation)
        raise SequenceParseError(f"Mode de lecture inconnu : {value!r} ({accepted})") from None


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
        char = token[index]
        if char == "[":
            end = token.index("]", index)  # présent : crochets déjà vérifiés
            residues.append(library.by_code(token[index + 1 : end]))
            index = end + 1
            continue
        if char.islower():
            raise SequenceParseError(
                f"Minuscule « {char} » en position {len(residues) + 1} : en code à une lettre, "
                "une minuscule désigne d'ordinaire un acide aminé D (le f de c(RGDfK) est la "
                "D-phénylalanine), que la bibliothèque ne décrit pas encore ; écrire la "
                "séquence en majuscules pour des acides aminés L"
            )
        residues.append(library.by_one_letter(char))
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


def _read_one_letter(tokens: list[str], library: ResidueLibrary) -> list[Residue]:
    if len(tokens) > 1:
        raise SequenceParseError(
            "En code à une lettre, les résidus s'écrivent d'un seul bloc, sans tiret entre "
            f"eux (ex. Ac-RGDS-NH2) ; ici {len(tokens)} blocs : {', '.join(tokens)}. "
            "Si la séquence est en code à trois lettres, choisir ce mode de lecture"
        )
    token = tokens[0]
    try:
        return _one_letter_block(token, library)
    except (SequenceParseError, UnknownResidueError) as error:
        if "[" in token or not library.has_code(token):
            raise
        hint = (
            f"{error}. « {token} » est aussi un code à trois lettres ({token.upper()}) : "
            "pour ce résidu seul, choisir le code à trois lettres"
        )
        raise type(error)(hint) from None


def _read_three_letter(tokens: list[str], library: ResidueLibrary) -> list[Residue]:
    token = tokens[0]
    if len(tokens) > 1 or "[" in token:
        return [_three_letter_token(t, library) for t in tokens]
    try:
        return [library.by_code(token)]
    except UnknownResidueError as error:
        raise UnknownResidueError(
            f"{error}. En code à trois lettres, chaque résidu occupe un bloc, séparé du "
            "suivant par un tiret (ex. Ac-Arg-Gly-Asp-Ser-NH2) ; si la séquence est en code "
            "à une lettre, choisir ce mode de lecture"
        ) from None


def parse_sequence(
    text: str,
    library: ResidueLibrary | None = None,
    *,
    notation: SequenceNotation | str,
) -> Peptide:
    """Construit un `Peptide` à partir d'une séquence texte, lue dans le code déclaré."""
    mode = _as_notation(notation)
    lib = library or load_default_library()
    cleaned = "".join(text.split())
    if not cleaned:
        raise SequenceParseError("Séquence vide")
    _check_brackets(cleaned)
    tokens = cleaned.split("-")
    if any(not t for t in tokens):
        raise SequenceParseError(f"Tiret mal placé dans la séquence : {text!r}")

    metadata: dict[str, str] = {"input": text, SEQUENCE_NOTATION: mode.value}

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

    if mode is SequenceNotation.ONE_LETTER:
        residues = _read_one_letter(tokens, lib)
    else:
        residues = _read_three_letter(tokens, lib)

    return Peptide(
        residues=tuple(residues),
        n_terminus=lib.n_terminus(n_alias),
        c_terminus=lib.c_terminus(c_alias),
        metadata=metadata,
    )


def _not_writable(peptide: Peptide) -> list[str]:
    """Déclarations qu'une séquence écrite ne sait pas porter."""
    lost = [link.describe() for link in peptide.links]
    lost += [modification.describe() for modification in peptide.modifications]
    lost += [
        f"état de protonation {residue.protonation} en position {position}"
        for position, residue in enumerate(peptide.residues, start=1)
        if residue.protonation is not None
    ]
    lost += [
        f"état de protonation {terminus.protonation} de l'extrémité {side}"
        for side, terminus in (("N", peptide.n_terminus), ("C", peptide.c_terminus))
        if terminus.protonation is not None
    ]
    return lost


def _residue_symbol(residue: Residue, mode: SequenceNotation) -> str:
    if not residue.standard or residue.one_letter is None:
        return f"[{residue.code}]"
    if mode is SequenceNotation.ONE_LETTER:
        return residue.one_letter
    return residue.code.title()


def format_sequence(
    peptide: Peptide,
    notation: SequenceNotation | str,
    library: ResidueLibrary | None = None,
) -> str:
    """Écrit `peptide` dans le code demandé, sous sa forme canonique.

    Lève `UnsupportedFeatureError` si le peptide porte ce qu'une séquence
    écrite ne représente pas (liaison hors squelette, modification, état de
    protonation), ou si le texte produit ne se relit pas comme le même peptide.
    """
    mode = _as_notation(notation)
    lib = library or load_default_library()
    if not peptide.residues:
        raise UnsupportedFeatureError("Peptide sans résidu : aucune séquence à écrire")
    lost = _not_writable(peptide)
    if lost:
        raise UnsupportedFeatureError(
            "Une séquence écrite ne porte ni liaison hors squelette, ni modification, ni "
            "état de protonation : écriture refusée plutôt que de perdre " + ", ".join(lost)
        )

    symbols = [_residue_symbol(residue, mode) for residue in peptide.residues]
    blocks: list[str] = []
    if N_TERMINUS_DEFAULTED not in peptide.metadata:
        blocks.append(lib.n_terminus_symbol(peptide.n_terminus.key))
    if mode is SequenceNotation.ONE_LETTER:
        blocks.append("".join(symbols))
    else:
        blocks.extend(symbols)
    if C_TERMINUS_DEFAULTED not in peptide.metadata:
        blocks.append(lib.c_terminus_symbol(peptide.c_terminus.key))
    text = "-".join(blocks)

    # Inventaire attendu contre texte produit : le texte est relu et comparé.
    try:
        reread = parse_sequence(text, lib, notation=mode)
    except PepForgeError as error:
        raise UnsupportedFeatureError(
            f"Écriture refusée : {text!r} ne se relit pas en {mode.label} ({error})"
        ) from None
    if reread != peptide:
        raise UnsupportedFeatureError(
            f"Écriture refusée : {text!r} se relit comme un autre peptide "
            f"({_describe(reread)}) que celui d'origine ({_describe(peptide)})"
        )
    return text


def convert_sequence(
    text: str,
    source: SequenceNotation | str,
    target: SequenceNotation | str,
    library: ResidueLibrary | None = None,
) -> str:
    """Lit `text` dans le code `source` et le réécrit dans le code `target`.

    `convert_sequence("Ac-RGDS-NH2", "one_letter", "three_letter")` donne
    `"Ac-Arg-Gly-Asp-Ser-NH2"`. Avec `source` égal à `target`, la séquence est
    seulement remise sous sa forme canonique (`NH2-RGDS-CO-NH2` → `H-RGDS-NH2`).
    """
    lib = library or load_default_library()
    return format_sequence(parse_sequence(text, lib, notation=source), target, lib)


def _describe(peptide: Peptide) -> str:
    return f"{peptide.n_terminus.key} / {peptide.three_letter_sequence} / {peptide.c_terminus.key}"
