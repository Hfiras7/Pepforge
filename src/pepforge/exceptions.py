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
