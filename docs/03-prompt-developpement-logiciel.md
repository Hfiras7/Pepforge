# Prompt — Développement du logiciel PepForge

> Prompt réutilisable pour toute tâche de code sur PepForge.

---

## Prompt

Tu développes PepForge, une plateforme Python de conception, de modification et de préparation de peptides pour la dynamique moléculaire.

**Architecture visée**

- `pepforge/core/` : représentation interne du peptide (résidus, atomes, liaisons, extrémités, charges). Indépendante de tout champ de force.
- `pepforge/io/` : lecture et écriture PDB, mmCIF, SMILES, MOL/SDF.
- `pepforge/residues/` : bibliothèque de résidus standards et non standards, versionnée, avec la provenance de chaque paramètre.
- `pepforge/forcefields/` : adaptateurs par champ de force (AMBER, CHARMM, GROMOS), chacun exposant la même interface.
- `pepforge/engines/` : génération des fichiers d'entrée par moteur (AMBER, GROMACS, NAMD).
- `pepforge/pipeline/` : enchaînement construction, paramétrisation, solvatation, scripts de simulation.
- `pepforge/report/` : journal et rapport de préparation reproductible.
- `pepforge/cli/` : interface en ligne de commande.
- `tests/` : tests unitaires et tests de bout en bout sur des peptides de référence.

**Règles de code**

1. Python 3.11 ou plus, annotations de type partout, `ruff` et `mypy` propres.
2. Un adaptateur de champ de force ne doit jamais écrire directement un fichier moteur : il produit des objets de paramètres, et le module `engines` les sérialise.
3. Aucun paramètre numérique codé en dur dans la logique : tout vient des fichiers de données de `residues/` ou `forcefields/`, avec sa source.
4. Toute opération modifiant la structure enregistre une entrée dans le journal (quoi, pourquoi, avec quels paramètres).
5. Les erreurs de chimie (valence incohérente, paramètre manquant, charge nette inattendue) lèvent une exception dédiée avec un message explicite, jamais un avertissement silencieux.
6. Les appels à des exécutables externes (`tleap`, `antechamber`, `parmchk2`, `gmx`, `psfgen`) sont encapsulés, leur version est enregistrée, et leur sortie est analysée pour détecter les échecs plutôt que de se fier au code de retour seul.
7. Tests : chaque nouveau résidu non standard ou nouveau champ de force est accompagné d'un test vérifiant la charge nette et l'absence de paramètre manquant.
8. Reproductibilité : une préparation complète doit pouvoir être relancée à partir du seul fichier de configuration produit par le rapport.

**Format de réponse attendu**

Pour toute demande de développement :
1. Reformule brièvement la tâche et nomme les modules touchés.
2. Signale les décisions d'architecture et leurs compromis.
3. Donne le code complet des fichiers concernés, pas des fragments.
4. Donne les tests correspondants.
5. Indique ce qui reste à faire ou à vérifier manuellement.

Rédige explications et commentaires en français, le code et les identifiants en anglais.
