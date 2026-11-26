"""Contrôles qualité exécutés à la sortie de chaque étape du pipeline.

Un contrôle est une fonction pure : il reçoit la sortie d'une étape et retourne
un `CheckResult`. Il ne modifie rien et n'appelle aucun outil externe. La
décision de bloquer le pipeline revient à `pepforge.pipeline`, qui lève
`QualityCheckError` dès qu'un résultat est en échec (`Status.FAIL`).

Les contrôles reprennent la liste de `docs/04-checklist-qualite.md` : chaque
contrôle porte l'identifiant de l'item qu'il couvre (`checklist_item`).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Generic, TypeVar

T_contra = TypeVar("T_contra", contravariant=True)


class Status(StrEnum):
    """Issue d'un contrôle."""

    PASS = "pass"
    """Contrôle satisfait."""

    WARN = "warn"
    """Point à relire par l'utilisateur, non bloquant (ex. extrémité par défaut)."""

    FAIL = "fail"
    """Contrôle en échec : le pipeline s'arrête."""


@dataclass(frozen=True)
class CheckResult:
    """Résultat d'un contrôle, conservé dans le rapport de préparation."""

    check_id: str
    status: Status
    message: str
    checklist_item: str | None = None
    """Référence à l'item de `docs/04-checklist-qualite.md` couvert, s'il y en a un."""

    details: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class Check(Generic[T_contra]):
    """Contrôle nommé, applicable à la sortie d'un type d'étape donné."""

    check_id: str
    description: str
    function: Callable[[T_contra], CheckResult]
    checklist_item: str | None = None

    def __call__(self, subject: T_contra) -> CheckResult:
        return self.function(subject)


def run_checks(checks: Iterable[Check[T_contra]], subject: T_contra) -> list[CheckResult]:
    """Exécute tous les contrôles, sans s'arrêter au premier échec."""
    return [check(subject) for check in checks]


def failures(results: Iterable[CheckResult]) -> list[CheckResult]:
    """Résultats en échec, dans l'ordre d'exécution."""
    return [result for result in results if result.status is Status.FAIL]
