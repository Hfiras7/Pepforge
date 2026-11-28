"""Exceptions dédiées de PepForge.

Règle du projet : une erreur de chimie ou de paramétrage lève toujours une
exception explicite, jamais un avertissement silencieux.
"""


class PepForgeError(Exception):
    """Classe de base de toutes les erreurs PepForge."""


class SequenceParseError(PepForgeError):
    """La séquence fournie ne peut pas être interprétée."""


class UnknownResidueError(PepForgeError):
    """Un résidu n'existe pas dans la bibliothèque chargée."""


class UnknownTerminusError(PepForgeError):
    """Une extrémité N- ou C-terminale n'est pas reconnue."""


class PeptideModelError(PepForgeError):
    """La représentation du peptide est incohérente.

    Position hors bornes, liaison d'un résidu à lui-même, liaison ou
    modification déclarée deux fois, clé interne mal formée, ou opération qui
    ferait disparaître en silence une liaison ou une modification déclarée.
    """


class UnsupportedFeatureError(PepForgeError):
    """Une déclaration du peptide n'est pas encore prise en charge par l'étape.

    Levée par une étape qui reçoit une déclaration (modification de chaîne
    latérale, liaison non linéaire) qu'elle ne sait pas traiter : le peptide
    est refusé plutôt que traité sans elle, ce qui serait une erreur silencieuse.
    """


class ChemistryError(PepForgeError):
    """Incohérence chimique (valence, stéréochimie, charge nette inattendue...)."""


class MissingParameterError(PepForgeError):
    """Un paramètre de champ de force est absent pour un résidu ou un fragment."""


class ReferenceDataError(PepForgeError):
    """Une donnée chimique de référence manque, est mal formée ou n'a pas de source.

    Par exemple un pKa à `null` dans `ionizable_groups.yaml`, un pKa dont la
    source n'est pas déclarée, un résidu non standard dont les groupes
    ionisables ne sont pas décrits, ou une extrémité que le fichier ne classe
    ni ionisable ni non ionisable. La donnée n'est jamais supposée : le fichier
    ou le peptide est refusé.
    """


class ConfigurationError(PepForgeError):
    """La configuration de préparation est incomplète ou incohérente."""


class ExternalToolError(PepForgeError):
    """Un exécutable externe (tleap, antechamber, gmx, psfgen...) a échoué."""


class StructureError(PepForgeError):
    """Structure 3D incohérente (atomes manquants, chevauchements, chaîne rompue...)."""


class StructureFormatError(StructureError):
    """Un fichier de structure est mal formé, d'un format non reconnu, ou lu de façon incomplète.

    Exemples : coordonnée illisible ou non finie, format annoncé par l'extension
    différent du contenu, nombre d'atomes lus différent du nombre d'atomes du
    fichier, modèle non choisi dans un fichier multi-modèles.
    """


class RetrievalError(StructureError):
    """Un fichier de structure n'a pas pu être obtenu (disque ou URL), ou il est refusé.

    Exemples : fichier absent, schéma d'URL refusé, délai dépassé, taille au-delà
    du plafond, réponse HTML, téléchargement tronqué.
    """


class QualityCheckError(PepForgeError):
    """Un contrôle qualité bloquant a échoué à la sortie d'une étape du pipeline."""

    def __init__(self, stage: str, failures: list[str]) -> None:
        self.stage = stage
        self.failures = failures
        details = "; ".join(failures)
        super().__init__(f"Étape « {stage} » : contrôle(s) en échec : {details}")


class ExtractionRefusedError(StructureError):
    """L'extraction d'un peptide depuis une chaîne de structure est refusée.

    Porte le rapport d'extraction complet (`report`), pour que le refus soit
    expliqué jonction par jonction et composant par composant.
    """

    def __init__(self, message: str, report: object) -> None:
        super().__init__(message)
        self.report = report
