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
- **Fichiers CHARMM** : un `.str` est un *stream file* portant topologie **et** paramètres, jamais des coordonnées. Les coordonnées sont dans le `.pdb` ou le `.crd`.
- **Sortie de CGenFF** : elle décrit une **molécule isolée**, pas un résidu de chaîne. Pour intégrer un résidu aza-β³, il faut en tirer un `RESI` ou un `PRES` avec ses atomes de liaison, puis ajuster les charges pour obtenir une charge entière compatible avec les résidus CHARMM36m voisins. Les pénalités CGenFF sont conservées et reportées.
- **Centre métallique** : ne jamais mélanger champ de force protéique CHARMM et paramètres métalliques AMBER dans un même système. Soit des paramètres compatibles CHARMM, soit tout le système marqué au format AMBER lu par NAMD, comme préparation distincte et vérifiée.

## Décisions encore ouvertes

Mise à jour du 26 novembre 2026, d'après le descriptif de thèse (voir [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), section 11) :

- **Moteur** : **NAMD** d'abord (affinement, dynamique moléculaire et minimisation ; VMD pour l'outillage). **GROMACS ensuite**, avec le même champ de force, comme conversion vérifiée. **AMBER en dernier** : un autre champ de force impose de dériver toute la bibliothèque aza-β³ une seconde fois.
- **Portée** : tout peptide de **2 à 100 résidus**, avec **une, deux ou plus** substitutions aza-β³, du même type ou de types différents. Pas de liste figée de séquences candidates.
- **Paramétrisation** : voie principale par **bibliothèque de résidus** (les résidus standards gardent CHARMM36m et sa correction CMAP) ; la voie « molécule entière dans CGenFF » sert de **contrôle croisé** sur un peptide court, pas de voie principale.
- **Système simulé** : **les deux**, le complexe peptide–récepteur étant le cas principal (docking puis dynamique).
- **Variante pseudopeptidique** : **aza-β³**, squelette `-NH-N(R)-CH2-CO-`. C'est la fonction centrale de l'outil (section 10 de l'architecture).
- **Champ de force** : **CHARMM36m + CGenFF**, confirmé par Firas. Modèle d'eau TIP3P modifié CHARMM. Fichiers à produire : `.rtf` (topologie de résidu), `.prm` (paramètres par type d'atome), `.str` (*stream file* : topologie **et** paramètres d'un résidu), `.psf` (système assemblé, via `psfgen`), coordonnées `.pdb`, et le script NAMD.
- **Radiomarquage** : les **deux** cas sont au périmètre, le peptide avec son chélateur sans métal **et** le peptide marqué, **quel que soit le métal**. Le mécanisme est générique (une entrée de données par métal), mais les **paramètres ne le sont jamais** : sans entrée validée et sourcée, `MissingParameterError`.
- **Intégration dans CHARMM36m** : c'est le vrai risque technique. Un aza-β³ modifie le **squelette**, pas seulement la chaîne latérale. Cinq familles de termes à vérifier (liaisons de jonction, angles, dièdres, impropres, non liées), **cinq patches `PRES` par type de résidu** (`standard → aza`, `aza → standard`, `aza → aza`, `N-ter aza`, `C-ter aza`), et énumération des termes depuis le `.psf` réellement produit. Voir [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), section 13.
- **CMAP** : le résidu aza-β³ n'en reçoit **aucun** (sa suite d'atomes de squelette n'est pas celle d'un acide aminé α). Les voisins gardent le leur, mais leur CMAP **traverse la jonction** : vérifier que ses références se résolvent, et qu'aucun terme n'est omis silencieusement. L'absence de CMAP est inscrite dans le rapport de chaque préparation concernée.
- **Métaux** : une entrée de données par couple métal–chélateur (`source` obligatoire), plus l'ajout manuel soumis aux mêmes contrôles. Ajouter un métal ne demande aucun code ; sans entrée validée, `MissingParameterError`.
- **Calcul quantique** : ORCA recommandé en premier, Psi4 en second (interface Python, automatisable), Gaussian seulement si licence déjà disponible. Les deux premiers derrière la même interface `QuantumTool`. Le **niveau de théorie est fixé une fois pour toute la bibliothèque** : deux entrées dérivées à des niveaux différents ne sont pas comparables.
- **Faisabilité sur poste de bureau** (section 16) : les balayages 1D de toute la bibliothèque tiennent en quelques jours ; une grille CMAP 2D sur **un** résidu modèle, transférée aux autres comme CHARMM le fait, tient en 1 à 2 jours ; une grille par entrée demanderait des mois. Le squelette aza-β³ a **trois** torsions, donc l'équivalent exact de CMAP serait 3D, ce que le format CMAP (grille 2D sur huit atomes) ne permet pas : vérifier ce point avant toute dépense de calcul. Le gros poste de calcul du projet n'est pas le quantique mais la dynamique, où **un GPU compte bien plus que les cœurs**.
- **Restent ouverts** : le logiciel quantique et le niveau de théorie à arrêter avec Firas, l'approche de paramétrisation du centre métallique, et le périmètre des balayages (unidimensionnels par entrée, approfondis sur un ou deux résidus de référence).

Malgré ces choix, le code garde moteur, champ de force et objectif **configurables** (`pepforge.config.PrepConfig`) et n'en suppose jamais un par défaut. Si une tâche dépend d'une décision encore ouverte, poser la question.

## Architecture du code

| Module | Rôle |
|---|---|
| `src/pepforge/core/` | Représentation interne (`Peptide`, `Structure`), indépendante de tout champ de force |
| `src/pepforge/config.py` | `PrepConfig`, `SolvationSpec` : choix de l'utilisateur, sans valeur par défaut |
| `src/pepforge/tools/` | `ExternalTool` : exécutables externes (module feuille) |
| `src/pepforge/checks/` | Contrôles qualité appliqués à la sortie de chaque étape |
| `src/pepforge/io/` | Lecture et écriture PDB, mmCIF, SMILES, MOL/SDF, et analyse des séquences |
| `src/pepforge/residues/` | Bibliothèque de résidus versionnée, avec la provenance de chaque donnée (`residues/data/`) |
| `src/pepforge/build/` | Construction et correction 3D, conformes à `StructureBuilder` / `StructureModifier` |
| `src/pepforge/forcefields/` | Adaptateurs de champ de force, tous conformes à `ForceFieldAdapter` |
| `src/pepforge/engines/` | Écriture des fichiers d'entrée par moteur, conformes à `EngineWriter` |
| `src/pepforge/pipeline/` | Enchaînement des étapes (`Stage`, `run_stage`) et contrôles après chaque étape |
| `src/pepforge/report/` | Journal des modifications et rapport reproductible |
| `src/pepforge/cli/` | Interface en ligne de commande (`pepforge`) |
| `tests/` | Tests pytest ; `test_architecture.py` vérifie la matrice de dépendances, `contracts.py` fournit les tests de contrat |

L'architecture détaillée (contrats, matrice de dépendances, flux, stratégie de tests) est dans [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), et l'ordre de construction dans [docs/PLAN-CONSTRUCTION.md](docs/PLAN-CONSTRUCTION.md). Toute modification d'un contrat ou d'une dépendance entre modules met ces documents à jour.

## Règles de code

1. Python 3.11 ou plus, annotations de type partout, `ruff` et `mypy --strict` propres.
2. Un adaptateur de champ de force n'écrit jamais directement un fichier moteur : il produit des objets de paramètres, que le module `engines` sérialise.
3. Aucun paramètre numérique codé en dur dans la logique : tout vient des fichiers de données de `residues/` ou `forcefields/`, avec sa source.
4. Toute opération modifiant la structure enregistre une entrée dans le journal (`pepforge.report.Journal`) : quoi, pourquoi, avec quels paramètres.
5. Les erreurs de chimie (valence incohérente, paramètre manquant, charge nette inattendue) lèvent une exception dédiée de `pepforge.exceptions`, jamais un avertissement silencieux.
6. Les appels à `tleap`, `antechamber`, `parmchk2`, `gmx`, `psfgen` passent par `pepforge.tools.ExternalTool` : version enregistrée, sortie analysée pour détecter les échecs.
7. Chaque nouveau résidu non standard ou champ de force est accompagné d'un test vérifiant la charge nette et l'absence de paramètre manquant. Chaque implémentation d'interface passe sa fonction de contrat de `tests/contracts.py`.
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
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) : modules, contrats, dépendances, flux et stratégie de vérification
- [docs/PLAN-CONSTRUCTION.md](docs/PLAN-CONSTRUCTION.md) : étapes de construction, critères d'acceptation et tests
