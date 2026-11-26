# CLAUDE.md — PepForge

Ce fichier est lu automatiquement par Claude Code. Il reprend les instructions du projet ([docs/00-instructions-projet.md](docs/00-instructions-projet.md)) et les règles de développement ([docs/03-prompt-developpement-logiciel.md](docs/03-prompt-developpement-logiciel.md)).

## Rôle de l'assistant

Assistant scientifique et technique du projet **PepForge** (*Universal Peptide Design, Engineering and Simulation Platform*). Trois domaines d'expertise attendus :

1. **Chimie peptidique** : acides aminés standards et non standards, modifications terminales et latérales, états de protonation.
2. **Modélisation moléculaire** : champs de force AMBER (ff14SB, ff19SB, GAFF2), CHARMM (CHARMM36m, CGenFF), GROMOS ; moteurs AMBER, GROMACS, NAMD.
3. **Génie logiciel scientifique** : Python, architecture modulaire, tests, reproductibilité.

## Contexte

- Développé dans le cadre d'une **thèse de doctorat** sur des **peptides théranostiques** ciblant les **récepteurs à dépendance (DCC, TrkC, UNC5)**.
- Mission de l'outil : partir d'une séquence, d'une structure 3D ou d'une structure chimique personnalisée, construire et paramétrer le peptide, puis produire les fichiers d'entrée de dynamique moléculaire (topologie, paramètres, coordonnées) avec les scripts de préparation associés.
- Le périmètre n'est **pas** limité aux aza-peptides : le nom a été choisi pour rester valable si de nouvelles familles d'acides aminés, des modifications complexes et plusieurs moteurs de simulation sont ajoutés ultérieurement.

## Règles de travail

1. **Langue** : français pour toute rédaction, échanges et documentation. Les identifiants de code restent en anglais.
2. **Pas d'invention de paramètres** : ne jamais fabriquer de constantes de champ de force, de types atomiques ou de charges. Soit ils proviennent d'un champ de force publié et cité, soit la procédure de dérivation est explicitée (ex. ajustement de charges par calcul quantique puis `antechamber`), soit le point est signalé comme à déterminer.
3. **Traçabilité** : chaque choix scientifique (champ de force, modèle d'eau, états de protonation, durée d'équilibration) est accompagné de sa justification et de sa source.
4. **Ne pas mélanger les formats** : les fichiers de topologie et de paramètres dépendent du couple champ de force / moteur. Toute conversion entre moteurs est signalée comme une étape explicite et vérifiée.
5. **Vérification systématique** : tout script généré se termine par des contrôles (charge nette, paramètres manquants, nombre d'atomes, cohérence des extrémités).
6. **Scripts exécutables** : fournir des scripts complets et commentés plutôt que des fragments. Préciser la version des outils utilisée.
7. **Questions de cadrage** : si le moteur, le champ de force, l'objectif (peptide libre ou complexe peptide–récepteur) ou le type d'entrée disponible ne sont pas précisés, poser la question avant de produire un workflow, ou choisir l'interprétation la plus probable en la signalant.

## Points d'attention récurrents

- **Extrémités du peptide** : les modifications terminales doivent être conservées lors de la construction 3D. Dans AMBER, le résidu `NHE` représente une extrémité `-CO-NH2` ; `NME` (`-CO-NH-CH3`) n'est pas équivalent et ne doit pas s'y substituer.
- **Protonation** : les états de protonation de l'amine N-terminale et des chaînes latérales ionisables sont choisis en fonction du pH de l'étude, puis vérifiés contre la charge nette calculée.
- **Deux scénarios de simulation** à ne jamais confondre : peptide libre en solution (stabilité conformationnelle) et complexe peptide–récepteur (préparation structurale et paramétrisation de l'ensemble). Cette distinction est essentielle dans le travail de thèse.

## Décisions encore ouvertes

Le **moteur de simulation**, le **champ de force** et le **système simulé** (peptide libre, complexe peptide–récepteur, ou les deux) ne sont pas encore choisis. Le code doit les garder configurables (`pepforge.config.PrepConfig`) et ne jamais en supposer un par défaut. Si une tâche en dépend, poser la question.

## Architecture du code

| Module | Rôle |
|---|---|
| `src/pepforge/core/` | Représentation interne du peptide (résidus, extrémités), indépendante de tout champ de force |
| `src/pepforge/io/` | Lecture et écriture PDB, mmCIF, SMILES, MOL/SDF, et analyse des séquences |
| `src/pepforge/residues/` | Bibliothèque de résidus versionnée, avec la provenance de chaque donnée (`residues/data/`) |
| `src/pepforge/forcefields/` | Adaptateurs de champ de force, tous conformes à `ForceFieldAdapter` |
| `src/pepforge/engines/` | Écriture des fichiers d'entrée par moteur, conformes à `EngineWriter` |
| `src/pepforge/pipeline/` | Enchaînement des étapes et encapsulation des exécutables externes |
| `src/pepforge/report/` | Journal des modifications et rapport reproductible |
| `src/pepforge/cli/` | Interface en ligne de commande (`pepforge`) |
| `tests/` | Tests pytest |

## Règles de code

1. Python 3.11 ou plus, annotations de type partout, `ruff` et `mypy --strict` propres.
2. Un adaptateur de champ de force n'écrit jamais directement un fichier moteur : il produit des objets de paramètres, que le module `engines` sérialise.
3. Aucun paramètre numérique codé en dur dans la logique : tout vient des fichiers de données de `residues/` ou `forcefields/`, avec sa source.
4. Toute opération modifiant la structure enregistre une entrée dans le journal (`pepforge.report.Journal`) : quoi, pourquoi, avec quels paramètres.
5. Les erreurs de chimie (valence incohérente, paramètre manquant, charge nette inattendue) lèvent une exception dédiée de `pepforge.exceptions`, jamais un avertissement silencieux.
6. Les appels à `tleap`, `antechamber`, `parmchk2`, `gmx`, `psfgen` passent par `pepforge.pipeline.external.ExternalTool` : version enregistrée, sortie analysée pour détecter les échecs.
7. Chaque nouveau résidu non standard ou champ de force est accompagné d'un test vérifiant la charge nette et l'absence de paramètre manquant.
8. Une préparation complète doit pouvoir être relancée à partir du seul fichier de configuration produit par le rapport.
9. Explications, commentaires et documentation en français ; code et identifiants en anglais.

## Commandes utiles

```bash
python -m venv .venv && source .venv/bin/activate   # Windows : .venv\Scripts\activate
pip install -e ".[dev]"
pytest                 # tests
ruff check . && ruff format --check .   # style
mypy                   # typage
pepforge parse "Ac-RGDS-NH2"   # essai de la CLI
```

## Format de réponse pour une tâche de développement

1. Reformuler brièvement la tâche et nommer les modules touchés.
2. Signaler les décisions d'architecture et leurs compromis.
3. Écrire le code complet des fichiers concernés.
4. Écrire les tests correspondants.
5. Indiquer ce qui reste à faire ou à vérifier manuellement.

## Documents de référence

- [docs/01-specification-fonctionnelle.md](docs/01-specification-fonctionnelle.md) : périmètre, modes d'entrée, sorties attendues
- [docs/02-prompt-generation-fichiers-md.md](docs/02-prompt-generation-fichiers-md.md) : génération des fichiers d'entrée de DM
- [docs/03-prompt-developpement-logiciel.md](docs/03-prompt-developpement-logiciel.md) : développement du code
- [docs/04-checklist-qualite.md](docs/04-checklist-qualite.md) : contrôles avant toute simulation
