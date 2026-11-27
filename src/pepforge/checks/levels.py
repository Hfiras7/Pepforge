"""Niveaux de vérification demandés, et audit de ce qui a réellement été exécuté.

Ce module répond à une exigence précise : l'utilisateur choisit jusqu'où pousser
la vérification d'une sortie (`ValidationLevel`), et **à chaque contrôle du
niveau demandé qui n'a pas pu être exécuté, l'outil signale une erreur**. C'est
le pendant, à l'échelle du protocole entier, de la règle du module
`pepforge.checks.inventory` : « je n'ai pas pu contrôler » ne se lit jamais
comme « tout va bien ».

Sans ce module, un niveau de validation serait une promesse : l'utilisateur
demanderait le niveau 3, l'outil exécuterait ce qu'il sait faire, et le rapport
ne dirait pas que les profils de torsion n'ont jamais été calculés. Ici, chaque
contrôle du catalogue de `docs/ARCHITECTURE.md` section 18 est déclaré avec son
niveau et son état de disponibilité, et `audit()` compare le catalogue aux
résultats réellement produits.

Le module est pur : il ne lance aucun contrôle et n'appelle aucun outil externe.
Il ne fait que confronter ce qui était demandé à ce qui a été rendu.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from pepforge.checks.base import CheckResult, Status
from pepforge.config import ValidationLevel


@dataclass(frozen=True)
class Control:
    """Un contrôle du catalogue de la section 18, avec son niveau et son état."""

    control_id: str
    """Identifiant, identique au `check_id` du résultat qui le satisfait."""

    level: ValidationLevel
    label: str

    available: bool
    """`True` si PepForge sait exécuter ce contrôle aujourd'hui."""

    pending_step: str | None = None
    """Étape de `docs/PLAN-CONSTRUCTION.md` qui le rendra disponible, le cas échéant."""

    checklist_item: str | None = None


CONTROLS: tuple[Control, ...] = (
    # Niveau 0 : conformité des fichiers. Tous portés par `checks.inventory.compare`,
    # sauf la géométrie, qui demande l'étape de contrôles de structure.
    Control(
        "inventory.residue_sequence",
        ValidationLevel.FILE_CONFORMITY,
        "séquence et ordre des résidus",
        available=True,
        checklist_item="structure.residue_count",
    ),
    Control(
        "inventory.atom_counts",
        ValidationLevel.FILE_CONFORMITY,
        "nombre d'atomes par résidu",
        available=True,
        checklist_item="structure.missing_atoms",
    ),
    Control(
        "inventory.backbone_bonds",
        ValidationLevel.FILE_CONFORMITY,
        "liaisons du squelette, donc absence de rupture de chaîne",
        available=True,
    ),
    Control(
        "inventory.n_terminus",
        ValidationLevel.FILE_CONFORMITY,
        "extrémité N conforme au patch demandé",
        available=True,
    ),
    Control(
        "inventory.c_terminus",
        ValidationLevel.FILE_CONFORMITY,
        "extrémité C conforme au patch demandé",
        available=True,
    ),
    Control(
        "inventory.disulfides",
        ValidationLevel.FILE_CONFORMITY,
        "ponts disulfure effectivement présents",
        available=True,
        checklist_item="structure.links_and_modifications",
    ),
    Control(
        "inventory.net_charge",
        ValidationLevel.FILE_CONFORMITY,
        "charge nette relue dans le fichier produit",
        available=True,
        checklist_item="parametrisation.net_charge",
    ),
    Control(
        "inventory.cmap_terms",
        ValidationLevel.FILE_CONFORMITY,
        "nombre de termes CMAP, aucun omis silencieusement",
        available=True,
    ),
    Control(
        "structure.geometry",
        ValidationLevel.FILE_CONFORMITY,
        "aucune coordonnée aberrante ni recouvrement d'atomes",
        available=False,
        pending_step="4 (contrôles de structure)",
        checklist_item="structure.clashes",
    ),
    # Niveau 1 : cohérence numérique. Demande le moteur installé.
    Control(
        "energy.single_point",
        ValidationLevel.NUMERIC_COHERENCE,
        "énergie en point unique, chargement du système sans paramètre manquant",
        available=False,
        pending_step="7 (topologie pour NAMD)",
    ),
    Control(
        "energy.minimization",
        ValidationLevel.NUMERIC_COHERENCE,
        "minimisation courte qui converge",
        available=False,
        pending_step="9 (protocole de simulation)",
    ),
    Control(
        "geometry.after_minimization",
        ValidationLevel.NUMERIC_COHERENCE,
        "géométrie conservée après minimisation",
        available=False,
        pending_step="9 (protocole de simulation)",
    ),
    # Niveau 2 : cohérence structurale.
    Control(
        "dynamics.short_stability",
        ValidationLevel.STRUCTURAL_COHERENCE,
        "dynamique courte : liaison N–N stable, ω proche de la planéité",
        available=False,
        pending_step="9 (protocole de simulation)",
    ),
    Control(
        "geometry.vs_quantum_fragment",
        ValidationLevel.STRUCTURAL_COHERENCE,
        "géométrie du résidu aza comparée au fragment capé optimisé en quantique",
        available=False,
        pending_step="11D (dérivation quantique)",
    ),
    # Niveau 3 : cohérence énergétique.
    Control(
        "torsion.profile_vs_quantum",
        ValidationLevel.ENERGETIC_COHERENCE,
        "profils de torsion du champ de force contre quantique, dont N–N",
        available=False,
        pending_step="11D (dérivation quantique)",
    ),
    Control(
        "dynamics.replicas_converge",
        ValidationLevel.ENERGETIC_COHERENCE,
        "répliques indépendantes convergeant vers les mêmes populations",
        available=False,
        pending_step="14 (analyses)",
    ),
    Control(
        "dynamics.control_peptide",
        ValidationLevel.ENERGETIC_COHERENCE,
        "peptide témoin non modifié, en CHARMM36m pur, mêmes conditions",
        available=False,
        pending_step="14 (analyses)",
    ),
    # Niveau 4 : accord expérimental. Hors de l'outil par nature.
    Control(
        "experiment.nmr",
        ValidationLevel.EXPERIMENTAL_AGREEMENT,
        "confrontation aux données RMN du peptide",
        available=False,
        pending_step="hors outil : confrontation manuelle",
    ),
    Control(
        "experiment.affinity",
        ValidationLevel.EXPERIMENTAL_AGREEMENT,
        "confrontation aux mesures d'affinité et de spécificité du complexe",
        available=False,
        pending_step="hors outil : confrontation manuelle",
    ),
)


def requested(level: ValidationLevel) -> tuple[Control, ...]:
    """Contrôles exigés par un niveau, les niveaux inférieurs compris.

    Les niveaux sont cumulatifs : un niveau 3 qui n'aurait pas vérifié la
    conformité des fichiers ne voudrait rien dire.
    """
    return tuple(c for c in CONTROLS if c.level <= level)


def unavailable(level: ValidationLevel) -> tuple[Control, ...]:
    """Contrôles du niveau demandé que PepForge ne sait pas encore exécuter.

    Permet d'avertir l'utilisateur **avant** de lancer une préparation longue,
    plutôt que de le découvrir dans le rapport final.
    """
    return tuple(c for c in requested(level) if not c.available)


def audit(level: ValidationLevel, produced: Iterable[CheckResult]) -> list[CheckResult]:
    """Confronte les contrôles exigés par `level` aux résultats réellement produits.

    Retourne un résultat par contrôle exigé :

    - le résultat produit lui-même, s'il existe ;
    - un `FAIL` « pas encore implémenté » si l'outil ne sait pas le faire ;
    - un `FAIL` « demandé mais non exécuté » sinon, c'est-à-dire quand l'outil
      saurait le faire mais que l'étape ne l'a pas lancé.

    Les deux derniers cas sont le cœur du mécanisme : demander un niveau et ne
    pas l'obtenir est une erreur, pas un silence.
    """
    by_id = {result.check_id: result for result in produced}
    results: list[CheckResult] = []
    for control in requested(level):
        existing = by_id.get(control.control_id)
        if existing is not None:
            results.append(existing)
            continue
        if not control.available:
            reason = (
                f"pas encore implémenté dans PepForge ({control.pending_step})"
                if control.pending_step
                else "pas encore implémenté dans PepForge"
            )
        else:
            reason = "demandé mais non exécuté par cette préparation"
        results.append(
            CheckResult(
                check_id=control.control_id,
                status=Status.FAIL,
                message=f"niveau {int(control.level)} — {control.label} : {reason}",
                checklist_item=control.checklist_item,
                details={
                    "niveau_demandé": str(int(level)),
                    "niveau_du_contrôle": str(int(control.level)),
                    "disponible": str(control.available),
                },
            )
        )
    return results
