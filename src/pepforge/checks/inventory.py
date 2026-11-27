"""Inventaire attendu contre fichier produit.

C'est le mécanisme central de la politique d'absence d'erreur silencieuse
(`docs/ARCHITECTURE.md`, section 17). Le principe : ne jamais croire l'outil
externe sur parole. On construit l'inventaire **attendu** depuis la
représentation interne, on relit le fichier **réellement produit**, et on
compare terme à terme.

Deux décisions de conception portent toute la valeur de ce module :

1. **Un champ non vérifiable est un échec, jamais un succès.** Si l'une des
   deux faces de la comparaison ignore une valeur, le résultat est `FAIL`
   avec la mention « non vérifiable », et non `PASS` ni un silence. C'est
   exactement ce qui distingue cette approche d'une vérification cosmétique :
   « je n'ai pas pu contrôler » ne doit jamais se lire comme « tout va bien ».
2. **La comparaison ne s'arrête pas au premier écart.** Elle retourne un
   résultat par champ, pour qu'un rapport montre tous les problèmes d'un coup
   plutôt que de les faire découvrir un par un.

Ce module ne lit aucun fichier : l'extraction de l'inventaire d'un `.psf`,
d'un `.prmtop` ou d'un `.top` est propre à chaque moteur et vit dans
`pepforge.engines`. Ici ne se trouve que la comparaison, qui est pure et donc
testable sans aucun outil externe installé.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pepforge.checks.base import CheckResult, Status
from pepforge.core.peptide import Peptide

CHARGE_TOLERANCE = 1e-3
"""Écart maximal toléré sur la charge nette, en charge élémentaire.

Une redistribution de charges de jonction approximative se trahit par un écart
de quelques centièmes (section 17.2) : la tolérance doit rester bien en dessous.
"""


@dataclass(frozen=True)
class Inventory:
    """Ce qu'un système contient, du point de vue vérifiable.

    Un champ à `None` signifie « inconnu », et non « vide » ou « sans
    importance ». La comparaison le traite comme un échec de vérifiabilité.
    """

    residue_sequence: tuple[str, ...] | None = None
    """Codes des résidus, dans l'ordre de la chaîne."""

    n_terminus: str | None = None
    """Clé interne de l'extrémité N (ex. `acetyl`), pas un nom de champ de force."""

    c_terminus: str | None = None
    atom_counts: dict[str, int] | None = None
    """Nombre d'atomes par résidu, indexé par `<position>:<code>`."""

    backbone_bonds: int | None = None
    """Nombre de liaisons peptidiques du squelette, ruptures de chaîne comprises."""

    disulfides: tuple[tuple[int, int], ...] | None = None
    """Paires de positions liées par un pont disulfure, normalisées et triées."""

    net_charge: float | None = None
    """Charge nette. Du côté produit, elle est **relue dans le fichier**."""

    cmap_terms: int | None = None
    """Nombre de termes CMAP. Un résidu aza-β³ n'en reçoit aucun (section 13.3)."""

    source: str = ""
    """D'où vient cet inventaire : la représentation interne, ou le fichier relu."""

    metadata: dict[str, str] = field(default_factory=dict)


def expected_from_peptide(peptide: Peptide) -> Inventory:
    """Construit la part de l'inventaire attendu déductible du seul `Peptide`.

    Les champs qui dépendent de la bibliothèque de résidus ou de l'adaptateur de
    champ de force (`atom_counts`, `net_charge`, `cmap_terms`) restent à `None` :
    ils seront renseignés par les étapes correspondantes du plan. Tant qu'ils le
    sont, la comparaison les signale comme non vérifiables, ce qui est le
    comportement voulu.
    """
    return Inventory(
        residue_sequence=tuple(r.code for r in peptide.residues),
        n_terminus=peptide.n_terminus.key,
        c_terminus=peptide.c_terminus.key,
        backbone_bonds=max(len(peptide.residues) - 1, 0),
        disulfides=(),
        source="représentation interne",
    )


def _compare_sequences(expected: tuple[str, ...], produced: tuple[str, ...]) -> tuple[Status, str]:
    if expected == produced:
        return Status.PASS, f"{len(expected)} résidus, dans le même ordre"
    if len(expected) != len(produced):
        missing = sorted(set(expected) - set(produced))
        detail = f" ; absents du fichier : {missing}" if missing else ""
        return (
            Status.FAIL,
            f"{len(expected)} résidus attendus, {len(produced)} dans le fichier{detail}",
        )
    divergences = [
        f"position {i + 1} : {e} attendu, {p} trouvé"
        for i, (e, p) in enumerate(zip(expected, produced, strict=True))
        if e != p
    ]
    return Status.FAIL, "; ".join(divergences)


def _compare_atom_counts(expected: dict[str, int], produced: dict[str, int]) -> tuple[Status, str]:
    problems: list[str] = []
    for key in sorted(set(expected) | set(produced)):
        exp = expected.get(key)
        prod = produced.get(key)
        if exp is None:
            problems.append(f"{key} : absent de l'attendu, {prod} atomes dans le fichier")
        elif prod is None:
            problems.append(f"{key} : {exp} atomes attendus, résidu absent du fichier")
        elif exp != prod:
            problems.append(f"{key} : {exp} atomes attendus, {prod} trouvés")
    if problems:
        return Status.FAIL, "; ".join(problems)
    return Status.PASS, f"{sum(expected.values())} atomes, répartition conforme"


def compare(expected: Inventory, produced: Inventory) -> list[CheckResult]:
    """Compare deux inventaires et retourne un résultat par champ.

    Tout champ inconnu d'un côté ou de l'autre produit un `FAIL` : une
    vérification impossible n'est jamais rapportée comme un succès.
    """
    results: list[CheckResult] = []

    def record(
        check_id: str,
        exp: object,
        prod: object,
        status: Status,
        message: str,
        checklist_item: str | None = None,
    ) -> None:
        results.append(
            CheckResult(
                check_id=check_id,
                status=status,
                message=message,
                checklist_item=checklist_item,
                details={
                    "attendu": repr(exp),
                    "produit": repr(prod),
                    "source_attendu": expected.source,
                    "source_produit": produced.source,
                },
            )
        )

    def unverifiable(check_id: str, exp: object, prod: object, label: str) -> None:
        side = "l'attendu" if exp is None else "le fichier produit"
        record(
            check_id,
            exp,
            prod,
            Status.FAIL,
            f"{label} : non vérifiable, valeur inconnue dans {side}",
        )

    pairs: list[tuple[str, object, object, str]] = [
        (
            "inventory.residue_sequence",
            expected.residue_sequence,
            produced.residue_sequence,
            "Séquence des résidus",
        ),
        ("inventory.n_terminus", expected.n_terminus, produced.n_terminus, "Extrémité N"),
        ("inventory.c_terminus", expected.c_terminus, produced.c_terminus, "Extrémité C"),
        (
            "inventory.atom_counts",
            expected.atom_counts,
            produced.atom_counts,
            "Nombre d'atomes par résidu",
        ),
        (
            "inventory.backbone_bonds",
            expected.backbone_bonds,
            produced.backbone_bonds,
            "Liaisons du squelette",
        ),
        ("inventory.disulfides", expected.disulfides, produced.disulfides, "Ponts disulfure"),
        ("inventory.net_charge", expected.net_charge, produced.net_charge, "Charge nette"),
        ("inventory.cmap_terms", expected.cmap_terms, produced.cmap_terms, "Termes CMAP"),
    ]

    for check_id, exp, prod, label in pairs:
        if exp is None or prod is None:
            unverifiable(check_id, exp, prod, label)
            continue

        if check_id == "inventory.residue_sequence":
            assert isinstance(exp, tuple) and isinstance(prod, tuple)
            status, message = _compare_sequences(exp, prod)
            record(check_id, exp, prod, status, message, "structure.residue_count")
        elif check_id == "inventory.atom_counts":
            assert isinstance(exp, dict) and isinstance(prod, dict)
            status, message = _compare_atom_counts(exp, prod)
            record(check_id, exp, prod, status, message, "structure.missing_atoms")
        elif check_id == "inventory.net_charge":
            assert isinstance(exp, float | int) and isinstance(prod, float | int)
            delta = abs(float(exp) - float(prod))
            if delta <= CHARGE_TOLERANCE:
                record(
                    check_id,
                    exp,
                    prod,
                    Status.PASS,
                    f"charge nette {prod:+g} conforme",
                    "parametrisation.net_charge",
                )
            else:
                record(
                    check_id,
                    exp,
                    prod,
                    Status.FAIL,
                    f"charge nette {prod:+g} dans le fichier, {float(exp):+g} attendue "
                    f"(écart {delta:.4g}, tolérance {CHARGE_TOLERANCE:g})",
                    "parametrisation.net_charge",
                )
        elif exp == prod:
            record(check_id, exp, prod, Status.PASS, f"{label} : conforme ({prod!r})")
        else:
            record(
                check_id,
                exp,
                prod,
                Status.FAIL,
                f"{label} : {exp!r} attendu, {prod!r} dans le fichier",
            )

    return results
