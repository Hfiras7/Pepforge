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


class ChemistryError(PepForgeError):
    """Incohérence chimique (valence, stéréochimie, charge nette inattendue...)."""


class MissingParameterError(PepForgeError):
    """Un paramètre de champ de force est absent pour un résidu ou un fragment."""


class ConfigurationError(PepForgeError):
    """La configuration de préparation est incomplète ou incohérente."""


class ExternalToolError(PepForgeError):
    """Un exécutable externe (tleap, antechamber, gmx, psfgen...) a échoué."""


class StructureError(PepForgeError):
    """Structure 3D incohérente (atomes manquants, chevauchements, chaîne rompue...)."""


class QualityCheckError(PepForgeError):
    """Un contrôle qualité bloquant a échoué à la sortie d'une étape du pipeline."""

    def __init__(self, stage: str, failures: list[str]) -> None:
        self.stage = stage
        self.failures = failures
        details = "; ".join(failures)
        super().__init__(f"Étape « {stage} » : contrôle(s) en échec : {details}")
