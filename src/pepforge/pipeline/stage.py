"""Contrat d'une étape du pipeline et exécution contrôlée.

Chaque étape est un composant remplaçable qui transforme une entrée typée en
une sortie typée (`Stage[InT, OutT]`). `run_stage` l'exécute toujours de la
même façon :

1. exécution de l'étape ;
2. exécution de tous les contrôles déclarés pour sa sortie ;
3. conservation des résultats dans le contexte (pour le rapport) ;
4. arrêt par `QualityCheckError` si un contrôle est en échec.

Une étape ne connaît ni l'étape précédente ni la suivante : seuls les types
échangés les relient. C'est ce qui permet de modifier ou remplacer un module
sans toucher aux autres.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, TypeVar, runtime_checkable

from pepforge.checks.base import Check, CheckResult, Status, failures, run_checks
from pepforge.config import PrepConfig
from pepforge.exceptions import QualityCheckError
from pepforge.report.journal import Journal

InT = TypeVar("InT", contravariant=True)
OutT = TypeVar("OutT", covariant=True)
In = TypeVar("In")
Out = TypeVar("Out")


@dataclass
class StageContext:
    """État partagé entre les étapes d'une préparation."""

    config: PrepConfig
    workdir: Path
    journal: Journal = field(default_factory=Journal)
    check_results: dict[str, list[CheckResult]] = field(default_factory=dict)
    """Nom d'étape → résultats de ses contrôles, dans l'ordre d'exécution."""


@runtime_checkable
class Stage(Protocol[InT, OutT]):
    """Étape du pipeline : une entrée typée, une sortie typée, aucun effet caché."""

    name: str

    def run(self, data: InT, context: StageContext) -> OutT:
        """Produit la sortie de l'étape. Toute modification de structure est journalisée."""
        ...


def run_stage(
    stage: Stage[In, Out],
    data: In,
    context: StageContext,
    checks: Sequence[Check[Out]] = (),
) -> Out:
    """Exécute une étape puis ses contrôles ; lève `QualityCheckError` en cas d'échec."""
    output = stage.run(data, context)
    results = run_checks(checks, output)
    context.check_results[stage.name] = results
    warnings = [r.check_id for r in results if r.status is Status.WARN]
    context.journal.record(
        "stage_completed",
        "fin d'étape et contrôles associés",
        stage=stage.name,
        checks=len(results),
        warnings=warnings,
    )
    failed = failures(results)
    if failed:
        raise QualityCheckError(stage.name, [f"{r.check_id} : {r.message}" for r in failed])
    return output
