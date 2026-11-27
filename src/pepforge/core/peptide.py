"""Modèle de données du peptide.

Cette représentation est indépendante de tout champ de force : aucun nom de
résidu, d'extrémité, de tautomère ou de patch propre à un champ de force n'y
figure (`tests/test_peptide.py` le vérifie). La traduction vers la nomenclature
d'un champ de force relève exclusivement de `pepforge.forcefields`.

Trois conventions tiennent tout le module.

1. **Immuabilité.** Un `Peptide` ne change jamais. Chaque opération (mutation,
   changement d'extrémité, état de protonation, liaison, modification) renvoie
   un nouveau peptide et laisse l'original intact. C'est l'appelant qui
   journalise l'opération (`pepforge.report.Journal`) : `core` ne dépend
   d'aucun autre module que `exceptions`.
2. **Positions numérotées à partir de 1**, de l'extrémité N vers l'extrémité C,
   comme dans la notation des séquences et dans les fichiers de structure.
3. **Clés internes neutres.** Extrémités, états de protonation et modifications
   sont désignés par des identifiants en minuscules (`acetyl`, `his_delta`...),
   jamais par un nom de champ de force. Une clé mal formée lève
   `PeptideModelError`.

Les liaisons non linéaires (`Link`) et les modifications de chaîne latérale
(`Modification`) sont **déclarées** ici ; leur prise en charge par les étapes
aval viendra plus tard. D'ici là, une étape qui ne sait pas traiter une
déclaration refuse le peptide (`UnsupportedFeatureError`) au lieu de l'ignorer.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field, replace
from enum import StrEnum

from pepforge.exceptions import PeptideModelError

N_TERMINUS_DEFAULTED = "n_terminus_source"
"""Clé de `Peptide.metadata` présente quand l'extrémité N n'a pas été déclarée.

Sa valeur dit d'où vient l'extrémité appliquée (« défaut de la bibliothèque »).
Les contrôles la lisent pour signaler une extrémité supposée plutôt que choisie.
"""

C_TERMINUS_DEFAULTED = "c_terminus_source"
"""Même rôle que `N_TERMINUS_DEFAULTED`, pour l'extrémité C."""

_INTERNAL_KEY = re.compile(r"[a-z][a-z0-9_]*")
_RESIDUE_CODE = re.compile(r"[A-Z0-9_]+")


def _check_internal_key(value: object, what: str) -> None:
    if not isinstance(value, str) or _INTERNAL_KEY.fullmatch(value) is None:
        raise PeptideModelError(
            f"{what} : clé interne mal formée {value!r} (attendu : lettres minuscules, "
            "chiffres et soulignés, en commençant par une lettre, comme « his_delta »)"
        )


def _check_position(value: object, what: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise PeptideModelError(
            f"{what} : position {value!r} invalide (entier à partir de 1 attendu)"
        )


@dataclass(frozen=True)
class Terminus:
    """Extrémité N- ou C-terminale (libre, acétylée, amidée...)."""

    key: str
    """Identifiant interne stable, par exemple `free_amine` ou `amide`."""

    label: str
    """Libellé lisible en français."""

    def __post_init__(self) -> None:
        _check_internal_key(self.key, "Extrémité")


@dataclass(frozen=True)
class Residue:
    """Résidu de la chaîne peptidique."""

    code: str
    """Code de la bibliothèque, en majuscules : `ARG` pour un standard, `AIB`...

    Lettres, chiffres et soulignés seulement : un tiret ou un crochet romprait
    la notation des séquences, et une casse libre ferait de `Arg` et `ARG` deux
    résidus différents.
    """

    one_letter: str | None = None
    """Code à une lettre s'il existe ; `None` pour la plupart des non standards."""

    standard: bool = True

    protonation: str | None = None
    """État de protonation choisi, en clé interne neutre ; `None` tant qu'aucun
    état n'est fixé.

    Exemples de clés : `his_delta`, `his_epsilon`, `his_protonated`. Le modèle
    ne vérifie que la forme de la clé. Le sens de l'état et sa compatibilité
    avec le résidu relèvent de l'étape de protonation (étape 2 du plan), qui
    connaît les groupes ionisables et la source de leurs pKa.
    """

    def __post_init__(self) -> None:
        if not isinstance(self.code, str) or _RESIDUE_CODE.fullmatch(self.code) is None:
            raise PeptideModelError(
                f"Code de résidu mal formé : {self.code!r} (lettres majuscules, chiffres "
                "et soulignés seulement)"
            )
        if self.one_letter is not None and (
            not isinstance(self.one_letter, str)
            or len(self.one_letter) != 1
            or not "A" <= self.one_letter <= "Z"
        ):
            raise PeptideModelError(
                f"Résidu {self.code} : code à une lettre mal formé {self.one_letter!r} "
                "(une lettre majuscule attendue)"
            )
        if self.protonation is not None:
            _check_internal_key(self.protonation, f"Protonation du résidu {self.code}")


class LinkKind(StrEnum):
    """Nature d'une liaison covalente qui sort de la chaîne linéaire."""

    DISULFIDE = "disulfide"
    """Pont disulfure entre les soufres des chaînes latérales de deux résidus."""

    HEAD_TO_TAIL = "head_to_tail"
    """Liaison amide entre l'amine α du premier résidu et le carboxyle α du dernier."""

    @property
    def label(self) -> str:
        """Libellé lisible en français."""
        return _LINK_LABELS[self]

    @property
    def involves_side_chains(self) -> bool:
        """Vrai si la liaison passe par les chaînes latérales des résidus liés.

        Une mutation remplace la chaîne latérale : elle est donc refusée à une
        position engagée dans une liaison de ce genre.
        """
        return self is LinkKind.DISULFIDE


_LINK_LABELS: dict[LinkKind, str] = {
    LinkKind.DISULFIDE: "pont disulfure",
    LinkKind.HEAD_TO_TAIL: "cyclisation tête-queue",
}


@dataclass(frozen=True)
class Link:
    """Liaison covalente entre deux résidus, hors liaisons du squelette linéaire.

    Les positions sont numérotées à partir de 1 et rangées dans l'ordre
    croissant à la construction : `Link(LinkKind.DISULFIDE, 7, 2)` et
    `Link(LinkKind.DISULFIDE, 2, 7)` désignent la même liaison.
    """

    kind: LinkKind
    first: int
    second: int

    def __post_init__(self) -> None:
        try:
            kind = LinkKind(self.kind)
        except ValueError:
            known = ", ".join(k.value for k in LinkKind)
            raise PeptideModelError(
                f"Type de liaison inconnu : {self.kind!r} (types connus : {known})"
            ) from None
        _check_position(self.first, kind.label)
        _check_position(self.second, kind.label)
        if self.first == self.second:
            raise PeptideModelError(
                f"{kind.label} : un résidu ne peut pas être lié à lui-même (position {self.first})"
            )
        object.__setattr__(self, "kind", kind)
        if self.first > self.second:
            first, second = self.second, self.first
            object.__setattr__(self, "first", first)
            object.__setattr__(self, "second", second)

    @property
    def positions(self) -> tuple[int, int]:
        return (self.first, self.second)

    def describe(self) -> str:
        """Description lisible, par exemple « pont disulfure 2–7 »."""
        return f"{self.kind.label} {self.first}–{self.second}"


@dataclass(frozen=True)
class Modification:
    """Modification déclarée de la chaîne latérale d'un résidu.

    `key` est une clé interne neutre (par exemple `phosphorylation`). Aucune
    bibliothèque de modifications n'existe encore : le modèle ne vérifie que la
    forme de la clé et la position, et les contrôles signalent toute
    modification déclarée comme non vérifiable.
    """

    key: str
    position: int

    def __post_init__(self) -> None:
        _check_internal_key(self.key, "Modification")
        _check_position(self.position, f"Modification {self.key}")

    def describe(self) -> str:
        return f"modification {self.key} en position {self.position}"


class _ReadOnlyMetadata(Mapping[str, str]):
    """Métadonnées en lecture seule, que l'on peut copier et sérialiser.

    `types.MappingProxyType` serait plus court, mais il ne se copie pas avec
    `copy.deepcopy`, ne se sérialise pas avec `pickle`, et fait échouer
    `dataclasses.asdict` sur le peptide entier.
    """

    __slots__ = ("_items",)

    def __init__(self, items: Mapping[str, str]) -> None:
        self._items = dict(items)

    def __getitem__(self, key: str) -> str:
        return self._items[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __repr__(self) -> str:
        return repr(self._items)


def _link_order(link: Link) -> tuple[int, int, str]:
    return (link.first, link.second, link.kind.value)


def _modification_order(modification: Modification) -> tuple[int, str]:
    return (modification.position, modification.key)


def _check_links(links: tuple[Link, ...], length: int) -> None:
    seen: set[Link] = set()
    disulfide_partner: dict[int, Link] = {}
    for link in links:
        if not isinstance(link, Link):
            raise PeptideModelError(f"Liaison attendue, reçu {link!r}")
        if link.second > length:
            raise PeptideModelError(
                f"{link.describe()} : position {link.second} hors du peptide "
                f"(positions valides de 1 à {length})"
            )
        if link in seen:
            raise PeptideModelError(f"{link.describe()} déclarée deux fois")
        seen.add(link)
        if link.kind is LinkKind.HEAD_TO_TAIL and link.positions != (1, length):
            raise PeptideModelError(
                f"{link.describe()} : une cyclisation tête-queue relie la position 1 "
                f"à la dernière position ({length})"
            )
        if link.kind is LinkKind.DISULFIDE:
            for position in link.positions:
                other = disulfide_partner.get(position)
                if other is not None:
                    raise PeptideModelError(
                        f"Position {position} engagée dans deux ponts disulfure "
                        f"({other.describe()} et {link.describe()})"
                    )
                disulfide_partner[position] = link


def _check_modifications(modifications: tuple[Modification, ...], length: int) -> None:
    seen: set[Modification] = set()
    for modification in modifications:
        if not isinstance(modification, Modification):
            raise PeptideModelError(f"Modification attendue, reçu {modification!r}")
        if modification.position > length:
            raise PeptideModelError(
                f"{modification.describe()} : position hors du peptide "
                f"(positions valides de 1 à {length})"
            )
        if modification in seen:
            raise PeptideModelError(f"{modification.describe()} déclarée deux fois")
        seen.add(modification)


def _without(metadata: Mapping[str, str], key: str) -> dict[str, str]:
    return {k: v for k, v in metadata.items() if k != key}


@dataclass(frozen=True)
class Peptide:
    """Chaîne de résidus, deux extrémités explicites, liaisons et modifications déclarées.

    Invariants vérifiés à la construction, `PeptideModelError` sinon :

    - `residues`, `links` et `modifications` sont des tuples ; liaisons et
      modifications sont rangées dans un ordre canonique, si bien que l'ordre
      de déclaration ne change pas l'égalité ;
    - toute position déclarée est comprise entre 1 et `len(peptide)` ;
    - une cyclisation tête-queue relie exactement la première et la dernière
      position ;
    - une position n'appartient qu'à un seul pont disulfure ;
    - aucune liaison ni modification n'est déclarée deux fois.

    `metadata` porte la provenance (texte saisi, extrémités supposées). Elle
    est en lecture seule et **exclue de l'égalité** : deux notations de la même
    molécule (`Ac-RGDS-NH2`, `Ac-Arg-Gly-Asp-Ser-NH2`) donnent le même peptide.
    Les opérations la conservent, sauf la marque d'extrémité supposée, retirée
    quand l'extrémité est fixée explicitement ; l'historique des opérations est
    dans le journal de l'appelant.
    """

    residues: tuple[Residue, ...]
    n_terminus: Terminus
    c_terminus: Terminus
    links: tuple[Link, ...] = ()
    modifications: tuple[Modification, ...] = ()
    metadata: Mapping[str, str] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        residues = tuple(self.residues)
        for position, residue in enumerate(residues, start=1):
            if not isinstance(residue, Residue):
                raise PeptideModelError(f"Position {position} : résidu attendu, reçu {residue!r}")
        for side, terminus in (("N", self.n_terminus), ("C", self.c_terminus)):
            if not isinstance(terminus, Terminus):
                raise PeptideModelError(f"Extrémité {side} : extrémité attendue, reçu {terminus!r}")
        links = tuple(self.links)
        modifications = tuple(self.modifications)
        _check_links(links, len(residues))
        _check_modifications(modifications, len(residues))
        object.__setattr__(self, "residues", residues)
        object.__setattr__(self, "links", tuple(sorted(links, key=_link_order)))
        object.__setattr__(
            self, "modifications", tuple(sorted(modifications, key=_modification_order))
        )
        object.__setattr__(self, "metadata", _ReadOnlyMetadata(self.metadata))

    def __len__(self) -> int:
        return len(self.residues)

    @property
    def one_letter_sequence(self) -> str:
        """Séquence à une lettre ; `X` pour un résidu sans code à une lettre."""
        return "".join(r.one_letter or "X" for r in self.residues)

    @property
    def three_letter_sequence(self) -> str:
        return "-".join(r.code for r in self.residues)

    def composition(self) -> dict[str, int]:
        """Nombre d'occurrences de chaque résidu (code de la bibliothèque)."""
        return dict(Counter(r.code for r in self.residues))

    @property
    def is_cyclic(self) -> bool:
        """Vrai si une cyclisation tête-queue est déclarée."""
        return any(link.kind is LinkKind.HEAD_TO_TAIL for link in self.links)

    @property
    def disulfides(self) -> tuple[tuple[int, int], ...]:
        """Paires de positions liées par un pont disulfure, dans l'ordre canonique."""
        return tuple(link.positions for link in self.links if link.kind is LinkKind.DISULFIDE)

    def _index(self, position: int) -> int:
        if (
            isinstance(position, bool)
            or not isinstance(position, int)
            or not 1 <= position <= len(self.residues)
        ):
            raise PeptideModelError(
                f"Position {position!r} hors du peptide "
                f"(positions valides de 1 à {len(self.residues)})"
            )
        return position - 1

    def _replace_residue(self, index: int, residue: Residue) -> Peptide:
        residues = (*self.residues[:index], residue, *self.residues[index + 1 :])
        return replace(self, residues=residues)

    def residue_at(self, position: int) -> Residue:
        """Résidu à une position numérotée à partir de 1."""
        return self.residues[self._index(position)]

    def mutate(self, position: int, residue: Residue) -> Peptide:
        """Renvoie un nouveau peptide où `residue` remplace le résidu de `position`.

        Le nouveau résidu arrive avec son propre état de protonation : celui de
        l'ancien ne lui est pas transmis, puisqu'il n'avait de sens que pour lui.

        La mutation est refusée (`PeptideModelError`) à une position qui porte
        une modification ou un pont disulfure : elle remplacerait la chaîne
        latérale dont dépend la déclaration, qui deviendrait fausse sans que
        rien ne le signale. Retirer d'abord la déclaration (`without_link`,
        `without_modification`). Une cyclisation tête-queue passe par le
        squelette : elle n'empêche pas la mutation du premier ou du dernier résidu.
        """
        index = self._index(position)
        if not isinstance(residue, Residue):
            raise PeptideModelError(f"Mutation : résidu attendu, reçu {residue!r}")
        blocking = [
            link.describe()
            for link in self.links
            if link.kind.involves_side_chains and position in link.positions
        ]
        blocking += [m.describe() for m in self.modifications if m.position == position]
        if blocking:
            raise PeptideModelError(
                f"Mutation refusée en position {position} "
                f"({self.residues[index].code} → {residue.code}) ; déclarations qui en "
                f"dépendent : {', '.join(blocking)}. Retirer d'abord ces déclarations."
            )
        return self._replace_residue(index, residue)

    def with_n_terminus(self, terminus: Terminus) -> Peptide:
        """Renvoie un peptide dont l'extrémité N est `terminus`, désormais déclarée."""
        return replace(
            self, n_terminus=terminus, metadata=_without(self.metadata, N_TERMINUS_DEFAULTED)
        )

    def with_c_terminus(self, terminus: Terminus) -> Peptide:
        """Renvoie un peptide dont l'extrémité C est `terminus`, désormais déclarée."""
        return replace(
            self, c_terminus=terminus, metadata=_without(self.metadata, C_TERMINUS_DEFAULTED)
        )

    def with_protonation(self, position: int, state: str | None) -> Peptide:
        """Fixe l'état de protonation d'un résidu, ou l'efface avec `None`."""
        index = self._index(position)
        return self._replace_residue(index, replace(self.residues[index], protonation=state))

    def with_link(self, link: Link) -> Peptide:
        """Renvoie un peptide qui déclare en plus la liaison `link`."""
        return replace(self, links=(*self.links, link))

    def without_link(self, link: Link) -> Peptide:
        """Renvoie un peptide sans la liaison `link`, qui doit être déclarée."""
        if link not in self.links:
            raise PeptideModelError(f"{link.describe()} n'est pas déclarée sur ce peptide")
        return replace(self, links=tuple(existing for existing in self.links if existing != link))

    def with_modification(self, modification: Modification) -> Peptide:
        """Renvoie un peptide qui déclare en plus la modification `modification`."""
        return replace(self, modifications=(*self.modifications, modification))

    def without_modification(self, modification: Modification) -> Peptide:
        """Renvoie un peptide sans `modification`, qui doit être déclarée."""
        if modification not in self.modifications:
            raise PeptideModelError(f"{modification.describe()} n'est pas déclarée sur ce peptide")
        return replace(
            self,
            modifications=tuple(m for m in self.modifications if m != modification),
        )
