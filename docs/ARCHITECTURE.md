# PepForge — Architecture

> Document de référence de l'architecture. Il décrit les modules, ce que chacun a le droit de connaître, les types de données qu'ils échangent, l'enchaînement du pipeline, et la façon dont chaque étape est revue, contrôlée et testée.
>
> Le plan de construction étape par étape est dans [PLAN-CONSTRUCTION.md](PLAN-CONSTRUCTION.md). Les règles de développement sont dans [`CLAUDE.md`](../CLAUDE.md).

## Sommaire

1. [Principes](#1-principes)
    - [1.1 Style d'architecture : monolithe modulaire](#11-style-darchitecture--monolithe-modulaire)
2. [Vue d'ensemble des modules](#2-vue-densemble-des-modules)
3. [Règles de dépendance](#3-règles-de-dépendance)
4. [Types de données échangés](#4-types-de-données-échangés)
5. [Contrat de chaque module](#5-contrat-de-chaque-module)
6. [Flux du pipeline](#6-flux-du-pipeline)
7. [Vérification : revue, contrôles et tests](#7-vérification--revue-contrôles-et-tests)
8. [Modifier un module isolément](#8-modifier-un-module-isolément)
9. [Décisions ouvertes et compromis](#9-décisions-ouvertes-et-compromis)
10. [Cas d'usage central : substitution aza-β³](#10-cas-dusage-central--substituer-un-acide-aminé-par-un-résidu-aza-β³)
11. [Contraintes venues du projet de thèse](#11-contraintes-venues-du-projet-de-thèse)
12. [Portée confirmée le 26 novembre 2026](#12-portée-confirmée-le-26-novembre-2026)
13. [Intégrer un résidu aza-β³ dans CHARMM36m](#13-intégrer-un-résidu-aza-β³-dans-charmm36m)
14. [Métaux : données livrées et ajout manuel](#14-métaux--données-livrées-et-ajout-manuel)
15. [Calcul quantique : recommandation et interface](#15-calcul-quantique--recommandation-et-interface)
    - [15.4 Ajouter concrètement ORCA et Psi4 au projet](#154-ajouter-concrètement-orca-et-psi4-au-projet)
    - [15.5 Automatiser davantage : les briques existantes](#155-automatiser-davantage--les-briques-existantes)
    - [15.6 Travailler avec Psi4 seul](#156-travailler-avec-psi4-seul--oui-et-cest-défendable)
16. [Faisabilité des calculs sur un poste de bureau](#16-faisabilité-des-calculs-sur-un-poste-de-bureau)
17. [Politique d'absence d'erreur silencieuse](#17-politique-dabsence-derreur-silencieuse)
18. [Comment vérifier et valider une sortie](#18-comment-vérifier-et-valider-une-sortie)
    - [18.8 Exemple complet de référence quantique](#188-exemple-complet-de-référence-quantique)
    - [18.9 Choisir le niveau de validation dans l'outil](#189-choisir-le-niveau-de-validation-dans-loutil)
19. [Distribution : un seul code, deux exécutables](#19-distribution--un-seul-code-deux-exécutables)
20. [Entrées, traitement par lots et sorties proposées](#20-entrées-traitement-par-lots-et-sorties-proposées)
21. [Interface utilisateur](#21-interface-utilisateur)

---

## 1. Principes

| # | Principe | Conséquence concrète |
|---|---|---|
| P1 | **Un module, une responsabilité** | Chaque dossier de `src/pepforge/` a un rôle unique (section 5). |
| P2 | **Les modules ne se parlent que par des types déclarés** | Les échanges passent par les types de la section 4 (`Peptide`, `Structure`, `ParameterSet`, `EngineFiles`...). Aucun dictionnaire « fourre-tout ». |
| P3 | **Les dépendances vont dans un seul sens** | Une matrice autorise ou interdit chaque import (section 3). Elle est vérifiée par `tests/test_architecture.py` : un import interdit fait échouer la CI. |
| P4 | **Les interfaces sont des `Protocol`** | Un adaptateur, un moteur, un constructeur ou une étape se remplace sans toucher au reste, tant qu'il respecte son `Protocol` et passe ses tests de contrat (`tests/contracts.py`). |
| P5 | **Chaque étape est contrôlée à sa sortie** | `pipeline.run_stage` exécute l'étape, puis ses contrôles (`checks`) ; un contrôle en échec arrête tout (`QualityCheckError`). |
| P6 | **Rien n'est supposé** | Moteur, champ de force, modèle d'eau, objectif, pH, solvatation, mode de lecture d'une séquence : aucune valeur par défaut (`PrepConfig`, `SolvationSpec`, `SequenceNotation`). Aucun paramètre numérique codé en dur : les données viennent de fichiers sourcés. |
| P7 | **Tout est tracé** | Chaque modification de structure est journalisée (`Journal`), chaque outil externe enregistre sa version, chaque contrôle laisse un résultat dans le rapport. |
| P8 | **Indépendance vis-à-vis du champ de force et du moteur** | `core`, `io`, `residues`, `build` et `checks` ne connaissent aucune nomenclature de champ de force (`NHE`, `HIE`...) ni aucun format moteur (`.prmtop`, `.top`, `.psf`). |

### 1.1 Style d'architecture : monolithe modulaire

Oui, c'est bien un **monolithe modulaire**, et le choix est délibéré. Un seul paquet Python, un seul processus, une seule base de code versionnée ensemble ; mais des frontières internes explicites, déclarées et **vérifiées à la compilation du test** plutôt que seulement documentées :

- les frontières sont les modules de `src/pepforge/`, chacun avec une responsabilité unique (P1) ;
- les contrats sont des `Protocol` typés et les types de données échangés de la section 4 (P2, P4) ;
- la direction des dépendances est une matrice exécutable, `ALLOWED` dans `tests/test_architecture.py`, qui fait échouer la CI sur un import interdit, un import relatif ou un cycle (P3).

C'est ce qui distingue un monolithe modulaire d'un monolithe ordinaire : dans un monolithe ordinaire, rien n'empêche `io` d'importer `forcefields` « juste pour cette fois », et la structure se dissout en un an. Ici, cet import ne passe pas la CI.

**Pourquoi pas des microservices ou des plugins installables séparément.** PepForge est un outil en ligne de commande, mono-utilisateur, qui tourne sur un poste de bureau (section 16). Il n'a ni montée en charge indépendante par composant, ni équipes séparées, ni déploiement continu : les trois seules raisons qui justifieraient le coût d'un découpage en services. Un découpage apporterait de la sérialisation, des versions à synchroniser et des pannes réseau, pour rien.

**Les points d'extension, qui sont la vraie demande.** Les adaptateurs de champ de force, les générateurs par moteur, les constructeurs 3D et les logiciels de calcul quantique sont des `Protocol` avec une implémentation par outil. Ajouter GROMACS, AMBER, Psi4 ou un nouveau chélateur, c'est ajouter un fichier qui respecte le `Protocol` et passe sa fonction de contrat de `tests/contracts.py` : aucune modification du reste. C'est l'indépendance que Firas a demandée, obtenue sans frontière de processus.

**Si cela devait changer un jour.** Les modules les plus susceptibles de sortir du monolithe sont `tools` (lancer un calcul quantique ou une dynamique sur une machine distante) et `engines`. Tous deux sont déjà derrière un `Protocol` et déjà des feuilles ou presque de la matrice de dépendances : les remplacer par un client distant resterait un changement local. C'est la propriété qu'il faut préserver, et elle ne demande pas de découper aujourd'hui.

---

## 2. Vue d'ensemble des modules

```mermaid
flowchart TB
    subgraph L5["Couche 5 — interface"]
        cli["cli<br/>commande pepforge"]
    end
    subgraph L4["Couche 4 — orchestration"]
        pipeline["pipeline<br/>étapes, run_stage, contexte"]
    end
    subgraph L3["Couche 3 — moteurs"]
        engines["engines<br/>EngineWriter → EngineFiles"]
    end
    subgraph L2["Couche 2 — traitements"]
        io["io<br/>séquence, PDB, mmCIF, SMILES, SDF"]
        build["build<br/>construction et correction 3D"]
        forcefields["forcefields<br/>ForceFieldAdapter → ParameterSet"]
    end
    subgraph L1["Couche 1 — services"]
        residues["residues<br/>bibliothèque sourcée"]
        checks["checks<br/>contrôles qualité"]
        report["report<br/>journal, rapport"]
    end
    subgraph L0["Couche 0 — fondations"]
        core["core<br/>Peptide, Structure"]
        config["config<br/>PrepConfig, SolvationSpec"]
        tools["tools<br/>ExternalTool"]
        exceptions["exceptions"]
    end

    cli --> pipeline
    pipeline --> engines & io & build & forcefields
    engines --> forcefields
    io & build & forcefields --> residues
    io & build & forcefields & engines --> report
    build & forcefields & engines --> tools
    report --> checks
    residues & checks & report --> core
    checks & report & forcefields & engines --> config
```

Les flèches indiquent « importe ». Toutes les couches peuvent importer `exceptions` (non dessiné pour la lisibilité). Une flèche ne remonte jamais vers une couche supérieure.

| Module | Rôle | État |
|---|---|---|
| `exceptions` | Hiérarchie d'exceptions dédiées (`PepForgeError` et dérivées) | en place |
| `core` | Représentation interne : `Peptide`, `Residue`, `Terminus`, `Structure`, `Atom` | `Peptide` en place, `Structure` minimale |
| `config` | `PrepConfig`, `SolvationSpec`, `Engine`, `SimulationGoal` ; YAML aller-retour | en place |
| `tools` | `ExternalTool` : exécutable externe, version, détection d'échec dans la sortie | en place |
| `residues` | Bibliothèque de résidus et d'extrémités chargée depuis `residues/data/*.yaml`, avec provenance | acides aminés standards et extrémités usuelles |
| `checks` | Contrôles qualité : fonctions pures `sortie d'étape → CheckResult` | socle + 2 contrôles sur `Peptide` |
| `report` | `Journal` des opérations ; rapport de préparation reproductible | journal en place, rapport à écrire |
| `io` | Lecture et écriture : séquences, puis PDB, mmCIF, SMILES, MOL/SDF | séquences en place |
| `build` | Construction 3D depuis un `Peptide`, corrections d'une `Structure` | interfaces seulement |
| `forcefields` | Adaptateurs de champ de force → `ParameterSet` | interface seulement |
| `engines` | Écriture des fichiers d'un moteur, solvatation | interface seulement |
| `pipeline` | `Stage`, `StageContext`, `run_stage`, étapes concrètes | socle + étape « séquence » |
| `cli` | Commande `pepforge` | `parse`, `config-template` |

---

## 3. Règles de dépendance

Matrice des imports internes autorisés. La source de vérité exécutable est le dictionnaire `ALLOWED` de [`tests/test_architecture.py`](../tests/test_architecture.py) : les deux doivent rester identiques.

| Module ↓ peut importer → | exceptions | core | config | tools | residues | checks | report | io | build | forcefields | engines | pipeline |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **exceptions** | — | | | | | | | | | | | |
| **core** | ✓ | — | | | | | | | | | | |
| **config** | ✓ | | — | | | | | | | | | |
| **tools** | ✓ | | | — | | | | | | | | |
| **residues** | ✓ | ✓ | | | — | | | | | | | |
| **checks** | ✓ | ✓ | ✓ | | | — | | | | | | |
| **report** | ✓ | ✓ | ✓ | | | ✓ | — | | | | | |
| **io** | ✓ | ✓ | | | ✓ | | ✓ | — | | | | |
| **build** | ✓ | ✓ | | ✓ | ✓ | | ✓ | | — | | | |
| **forcefields** | ✓ | ✓ | ✓ | ✓ | ✓ | | ✓ | | | — | | |
| **engines** | ✓ | ✓ | ✓ | ✓ | | | ✓ | | | ✓ | — | |
| **pipeline** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **cli** | ✓ | ✓ | ✓ | | | ✓ | ✓ | ✓ | | | | ✓ |

Points essentiels :

- **`io`, `build` et `forcefields` ne se connaissent pas.** Ils sont reliés uniquement par `pipeline`, à travers les types de `core`. On peut donc remplacer le constructeur 3D sans toucher aux champs de force, et inversement.
- **`forcefields` n'importe jamais `engines`.** Un adaptateur produit des paramètres, il n'écrit pas de fichier moteur (règle 2 de `CLAUDE.md`). `engines` importe `forcefields` uniquement pour le type `ParameterSet`.
- **`tools` est une feuille.** `ExternalTool` est utilisé par `build`, `forcefields` et `engines` ; s'il était resté dans `pipeline`, ces modules auraient dû importer leur propre orchestrateur, ce qui crée une dépendance circulaire. C'est pourquoi il a été déplacé de `pepforge.pipeline.external` vers `pepforge.tools.external`.
- **`cli` ne contient aucune logique scientifique.** Il lit les arguments, appelle `pipeline` (ou `io` pour les commandes d'inspection) et affiche.
- Les imports relatifs (`from .x import y`) sont interdits, pour que la matrice reste vérifiable.
- **Bibliothèques tierces confinées** (étape 3). gemmi, qui lit et écrit PDB et mmCIF, n'est importé que par `io/structure_formats.py` ; le reste de PepForge ne voit que `Structure`. Changer de bibliothèque ne touche donc qu'un fichier. Le test `test_third_party_library_is_confined` vérifie cette règle, et toute nouvelle bibliothèque tierce s'y déclare avec le fichier qui a le droit de l'importer.

Le test vérifie aussi que la matrice ne contient pas de cycle, et que tout nouveau module y est déclaré.

---

## 4. Types de données échangés

Ce sont les seuls objets qui traversent les frontières entre modules. Chaque type appartient à un module unique, qui seul le modifie.

| Type | Module propriétaire | Produit par | Consommé par | Invariants |
|---|---|---|---|---|
| `PrepConfig` | `config` | utilisateur (YAML), `cli` | `pipeline`, `checks`, `report` | Aucun champ décisionnel n'a de valeur par défaut. `missing_decisions()` liste ce qui manque. Aller-retour YAML sans perte. |
| `SolvationSpec` | `config` | `PrepConfig.solvation` | `engines` | Aucune valeur par défaut ; vérifié par l'étape de solvatation. |
| `Peptide` | `core` | `io` (séquence, structure) ; ses propres opérations (mutation, extrémités, protonation, liaisons, modifications), appelées notamment par `build` | `build`, `forcefields`, `checks` | **Immuable** : chaque opération renvoie un nouveau peptide et l'appelant la journalise. Positions numérotées à partir de 1. Extrémités N et C toujours explicites ; une extrémité par défaut est signalée dans `metadata` (`N_TERMINUS_DEFAULTED`, `C_TERMINUS_DEFAULTED`). Liaisons et modifications rangées dans un ordre canonique, positions dans les bornes, aucun doublon. `metadata` en lecture seule et exclue de l'égalité. Aucun nom propre à un champ de force (vérifié par un test). |
| `Residue`, `Terminus` | `core` | `residues` | tous | Immuables ; identifiants internes stables (`amide`, `acetyl`...) ; code de résidu en majuscules ; `Residue.protonation` et, depuis l'étape 2, `Terminus.protonation` facultatifs, en clé interne neutre (`his_epsilon`, `nterm_protonated`...). |
| `ProtonationResult`, `GroupAssignment`, `IonizationForm` | `core` | `residues.assign_protonation`, via `pipeline.ProtonationStage` | `checks`, `report`, étapes aval | *Ajouté à l'étape 2.* Peptide protoné, pH, un `GroupAssignment` par groupe ionisable (site, état, forme, forme majoritaire, charge, pKa, source, imposé ou non), charge nette **entière** des états retenus et charge moyenne au pH. Placé dans `core` et non dans `residues` pour que `checks` puisse le contrôler sans dépendre de `residues` (section 3, matrice inchangée). |
| `Link`, `LinkKind`, `Modification` | `core` | `Peptide.with_link`, `Peptide.with_modification` | `checks` | **Déclarés**, pas encore pris en charge en aval : une étape qui ne sait pas les traiter lève `UnsupportedFeatureError`. Positions d'une liaison rangées dans l'ordre croissant ; une cyclisation tête-queue relie la première et la dernière position ; une position n'appartient qu'à un seul pont disulfure. |
| `Structure`, `Atom`, `ResidueId`, `Connection` | `core` | `io` (PDB, mmCIF), `build` | `build`, `engines`, `checks` | Immuables ; coordonnées en Å, toujours finies ; nomenclature PDB, jamais celle d'un champ de force. Depuis l'étape 3 : code d'insertion, occupation et facteur B par atome (`None` s'ils ne viennent pas d'un fichier, jamais remplis par défaut), modèle lu, liaisons déclarées par le fichier (`SSBOND`, `LINK`, `struct_conn`) et résidus déclarés non observés (`REMARK 465`, `pdbx_unobs_or_zero_occ_residues`). Un résidu se désigne par chaîne, numéro et code d'insertion de l'auteur. |
| `RetrievedFile`, `DownloadLimits` | `io` | `io.retrieval` | `io.structure` | Contenu exact reçu (disque ou URL), avec origine, taille, SHA-256 des octets reçus et date. Bornes de téléchargement techniques : délai, taille maximale, `http` sur demande seulement. |
| `StructureInventory` | `io` | `io.structure_inventory` | `cli`, rapport ; contrôles de l'étape 4 | Description sans modification : chaînes, eaux, ions, composants non reconnus, hydrogènes, conformations alternatives, codes d'insertion, trous de numérotation, résidus non observés, liaisons déclarées, et liste des décisions à prendre. |
| `ResidueLibrary` | `residues` | `residues.load_default_library()` | `io`, `build`, `forcefields` | Chaque fichier chargé déclare sa source (`sources`). |
| `ParameterSet`, `ParameterSource` | `forcefields` | un `ForceFieldAdapter` | `engines`, `checks`, `report` | `missing` vide avant toute écriture moteur ; `sources` non vide. |
| `EngineFiles` | `engines` | un `EngineWriter` | `pipeline`, `checks`, `report` | Charge nette et nombre d'atomes **relus** dans les fichiers produits, pas recopiés de l'entrée ; versions des outils présentes. |
| `ToolResult` | `tools` | `ExternalTool.run()` | `build`, `forcefields`, `engines` | Version de l'outil renseignée ; une sortie contenant un motif d'échec lève `ExternalToolError`. |
| `JournalEntry`, `Journal` | `report` | tout composant qui modifie une structure, `run_stage` | `report` | Ajout seulement, jamais de suppression. |
| `CheckResult` | `checks` | `Check` | `pipeline`, `report` | `FAIL` arrête le pipeline ; `WARN` est reporté dans le rapport. |
| `StageContext` | `pipeline` | `pipeline` | les étapes | Contient `config`, `journal`, `workdir`, résultats des contrôles. |

```mermaid
classDiagram
    class Peptide {
        residues: tuple~Residue~
        n_terminus: Terminus
        c_terminus: Terminus
        links: tuple~Link~
        modifications: tuple~Modification~
        metadata: Mapping en lecture seule
    }
    class Link {
        kind: LinkKind
        first: int
        second: int
    }
    class Modification {
        key: str
        position: int
    }
    class Structure {
        atoms: tuple~Atom~
        peptide: Peptide optionnel
        source: str
    }
    class ParameterSet {
        forcefield: str
        residue_names: list~str~
        sources: dict
        parameter_files: list~ParameterSource~
        net_charge: float optionnel
        missing: list~str~
    }
    class EngineFiles {
        engine: Engine
        files: dict
        solvated: bool
        net_charge: float optionnel
        atom_count: int optionnel
        tool_versions: dict
    }
    Peptide --> Link : déclare
    Peptide --> Modification : déclare
    Structure --> Peptide : décrit
    ParameterSet ..> Peptide : paramètre
    EngineFiles ..> Structure : sérialise
    EngineFiles ..> ParameterSet : sérialise
```

---

## 5. Contrat de chaque module

Pour chaque module : ce qu'il reçoit, ce qu'il rend, les erreurs qu'il lève, et ce qu'il ne doit jamais faire.

### 5.1 `core` — représentation interne

- **Entrées / sorties** : aucune entrée externe ; fournit `Peptide`, `Residue`, `Terminus`, `Link`, `LinkKind`, `Modification`, `Structure`, `Atom`, et depuis l'étape 2 `ProtonationResult`, `GroupAssignment`, `IonizationForm` (données seulement : les pKa vivent dans `residues/data/`).
- **Erreurs** : `PeptideModelError` (position hors bornes, liaison d'un résidu à lui-même, cyclisation tête-queue qui ne relie pas les deux extrémités, position engagée dans deux ponts disulfure, liaison ou modification déclarée deux fois, clé interne mal formée, mutation qui rendrait fausse une déclaration).
- **Garanties** (étape 1) : `Peptide` immuable ; `mutate`, `with_n_terminus`, `with_c_terminus`, `with_protonation`, `with_link`, `without_link`, `with_modification` et `without_modification`, et depuis l'étape 2 `with_n_terminus_protonation` et `with_c_terminus_protonation` (qui, eux, conservent la marque d'extrémité supposée), renvoient un nouveau peptide et laissent l'original intact ; positions numérotées à partir de 1. La journalisation revient à l'appelant, puisque `core` ne dépend que d'`exceptions`. Une mutation est **refusée** à une position qui porte un pont disulfure ou une modification, parce qu'elle remplacerait la chaîne latérale dont dépend la déclaration ; la cyclisation tête-queue, qui passe par le squelette, ne l'empêche pas.
- **Ne fait jamais** : lire un fichier, connaître un champ de force ou un moteur (un test cherche leurs noms de résidus et de patches dans `src/pepforge/core/`), contenir une constante physique, juger de la chimie d'une déclaration, qui relève de `checks`.
- **Évolutions prévues** : plusieurs chaînes (complexe) ; prise en charge des liaisons et des modifications par les étapes aval.

### 5.2 `config` — configuration d'une préparation

- **Entrées** : dictionnaire ou fichier YAML. **Sorties** : `PrepConfig`.
- **Erreurs** : `ConfigurationError` (clé inconnue, valeur invalide, décision manquante via `require_complete()`).
- **Garantie** : le YAML écrit par le rapport suffit à relancer la préparation (règle 8).
- **Étape 2** : champ `protonation` (site → état imposé, `"N"`, `"C"` ou une position), vide par défaut, ce qui ne suppose rien. Il est indispensable pour une histidine neutre, dont le pKa ne départage pas les tautomères.

### 5.3 `tools` — exécutables externes

- **Entrées** : nom de l'exécutable, arguments, motifs d'échec. **Sorties** : `ToolResult`.
- **Erreurs** : `ExternalToolError` (exécutable absent, code de retour non nul, motif d'échec dans la sortie).
- **Garantie** : la version est toujours enregistrée ; le code de retour seul ne suffit pas (tleap retourne 0 même en cas de paramètre manquant).

### 5.4 `residues` — bibliothèque sourcée

- **Entrées** : fichiers `residues/data/*.yaml`, chacun avec un champ `source`. **Sorties** : `ResidueLibrary`.
- **Erreurs** : `UnknownResidueError`, `UnknownTerminusError` ; depuis l'étape 2, `ReferenceDataError` (pKa manquant ou sans source, résidu ou extrémité non classé, fichier incohérent), `ChemistryError` (état qui n'appartient pas au groupe, état sur un site non ionisable, pH égal au pKa ou tautomère d'histidine sans état imposé), `UnsupportedFeatureError` (peptide cyclique ou modifié), `ConfigurationError` (pH hors de 0 à 14).
- **Étape 2** : `ionizable_groups.yaml` (groupes, états, charges formelles, pKa et sources, fenêtre d'alerte justifiée) et `residues.protonation` : `load_ionizable_groups`, `assign_protonation(peptide, ph) -> ProtonationResult`, `theoretical_net_charge(peptide, ph) -> float`. Chaque résidu standard et chaque extrémité est classé ionisable ou non : rien n'est supposé neutre faute d'être cité.
- **Ne fait jamais** : contenir des paramètres de champ de force (ceux-ci vivent dans `forcefields/data/`).

### 5.5 `checks` — contrôles qualité

- **Entrées** : la sortie d'une étape (`Peptide`, `Structure`, `ParameterSet`, `EngineFiles`). **Sorties** : `CheckResult` (`PASS`, `WARN`, `FAIL`).
- **Règles** : un contrôle est une fonction pure, sans effet de bord ni appel externe ; il référence l'item de [04-checklist-qualite.md](04-checklist-qualite.md) qu'il couvre (`checklist_item`).
- **Ne fait jamais** : décider d'arrêter le pipeline (c'est `pipeline.run_stage` qui le fait) ni corriger ce qu'il constate.
- **Erreurs** : `UnsupportedFeatureError`, levée par `expected_from_peptide` quand l'inventaire attendu ne sait pas représenter une déclaration du peptide (cyclisation tête-queue, modification de chaîne latérale). Ce n'est pas un contrôle qui arrête le pipeline : c'est un inventaire qui refuse d'être construit faux (section 17.1).

### 5.6 `report` — journal et rapport

- **Entrées** : opérations (`Journal.record`), résultats de contrôles, `PrepConfig`, versions d'outils. **Sorties** : journal YAML, rapport de préparation, configuration rejouable.
- **Garantie** : un rapport contient tout ce qu'il faut pour refaire la préparation à l'identique.

### 5.7 `io` — lecture et écriture

- **Entrées** : texte (séquence), fichiers PDB, mmCIF, SMILES, MOL/SDF. **Sorties** : `Peptide`, `Structure` (et, pour le mode « structure personnalisée », un fragment chimique à définir à l'étape 11 du plan).
- **Erreurs** : `SequenceParseError` (dont les crochets mal formés, le mode de lecture inconnu et la minuscule en code à une lettre), `UnknownResidueError`, `UnknownTerminusError`, `UnsupportedFeatureError` (écriture d'une séquence qui perdrait une déclaration ou ne se relirait pas à l'identique ; structure non représentable dans le format demandé), `StructureError`, et ses deux sous-classes de l'étape 3 : `StructureFormatError` (fichier mal formé, format contredit, lecture incomplète, modèle non choisi) et `RetrievalError` (fichier ou URL refusés).
- **Mode de lecture déclaré au départ** (revue du module 1) : `parse_sequence(texte, notation=...)`, avec `SequenceNotation.ONE_LETTER` (code à une lettre, monolitère) ou `SequenceNotation.THREE_LETTER` (code à trois lettres, trilitère). Le paramètre n'a **pas de valeur par défaut**, ni dans l'API, ni dans la CLI (`--notation 1` ou `3`), ni dans `SequenceInputStage`, et `PrepConfig.sequence_notation` est une décision exigée par `missing_decisions()` dès qu'une séquence est donnée. Rien n'est deviné : `ALA` est le tripeptide Ala-Leu-Ala en code à une lettre, l'alanine seule en code à trois lettres. Le mode est inscrit dans `Peptide.metadata` (`SEQUENCE_NOTATION`) et au journal. Une minuscule en code à une lettre est refusée, parce qu'elle désigne d'ordinaire un acide aminé D (le `f` de `c(RGDfK)`), que la bibliothèque ne décrit pas encore. Ce choix explicite remplace la lecture devinée signalée par avertissement, introduite puis retirée pendant la même revue.
- **Conversion entre les deux codes** (revue du module 1) : `convert_sequence(texte, source, cible)` lit la séquence dans son code et la réécrit dans l'autre, à la manière du réglage `seq_view_format` de PyMOL (0 pour les codes à une lettre, 1 pour les noms à trois lettres, [PyMOL Wiki](https://wiki.pymol.org/index.php?title=Seq_view_format)). `format_sequence(peptide, code)` produit la forme canonique : résidus standards selon l'IUPAC (`R`, `Arg`), non standards entre crochets avec leur code de bibliothèque (`[AIB]`), extrémités sous leur forme usuelle (premier alias de `residues/data/termini.yaml`) ; une extrémité non déclarée le reste. Le texte produit est **relu et comparé** au peptide (inventaire attendu contre sortie, section 17.1), et toute différence lève `UnsupportedFeatureError` : une conversion qui changerait la molécule n'est jamais rendue. Une séquence écrite ne porte ni liaison hors squelette, ni modification, ni état de protonation : leur présence est refusée, pas omise.
- **Notation des non standards** (étape 1) : code de la bibliothèque entre crochets, dans les deux notations (`Ac-RG[Aib]S-NH2`, `Ac-Arg-Gly-[Aib]-Ser-NH2`), convention reprise de HELM (Zhang *et al.*, *J. Chem. Inf. Model.* 2012, 52, 2796–2806). La règle de lecture complète est dans la docstring de [`io/sequence.py`](../src/pepforge/io/sequence.py).
- **Structures (étape 3, module 3)** : deux modes d'importation, une seule lecture (demande de Firas, 27 novembre 2026). `read_structure(chemin, journal=..., model=..., fmt=...)` lit un fichier local ; `read_structure_url(url, journal=..., limits=...)` le télécharge sous contrôle (`https` seulement sauf `allow_http`, aucun identifiant dans l'URL, délai, taille plafonnée pendant la lecture, réponse autre que 200 ou HTML refusée, téléchargement tronqué refusé, gzip décompressé sous le même plafond), puis le lit de la même façon. Le journal reçoit la provenance (chemin ou URL, URL finale, date, taille, SHA-256), le format, le modèle et la version de gemmi.
- **Registre de formats** (`io/structure_formats.py`) : chaque format inscrit son lecteur et son écrivain (`READERS`, `WRITERS`), ses extensions et sa signature de contenu. L'extension, le contenu et le format demandé doivent concorder, sinon `StructureFormatError`. Ajouter un format ne touche que ce fichier.
- **Lecture stricte, inventaire attendu contre lecture produite** : PepForge relit lui-même chaque enregistrement d'atome du fichier brut (colonnes du PDB 3.3, table `_atom_site` du mmCIF) et compare la structure de gemmi à cet inventaire, atome par atome (identité et coordonnées). Raison : gemmi 0.7.5 lit `1x.104` comme `1,0` en PDB, `?` comme NaN en mmCIF, et remplit une occupation absente par 1,0 et un facteur B absent par 20,0 (constaté le 27 novembre 2026). Atome en double, modèle sans `ENDMDL`, `NUMMDL` contredit par le nombre de modèles : refus.
- **Modèles** : un fichier multi-modèles exige `model=` (`--model` en ligne de commande) ; le premier n'est jamais pris en silence.
- **Écriture** (`write_structure`) : le fichier écrit sur le disque est relu et comparé atome par atome, liaisons déclarées comprises, avant d'être placé à son chemin final. Ce que le PDB ne sait pas porter (plus de 99 999 atomes, chaîne de deux caractères, numéro hors de [−999, 9 999], liaison hydrogène déclarée) lève `UnsupportedFeatureError` au lieu d'être tronqué.
- **Inventaire** (`io/structure_inventory.py`) : `inventory(structure)` décrit sans modifier et liste les décisions à prendre. Les eaux (`HOH`, `DOD`) et les groupes terminaux (`ACE`, `NH2`, `NME`) viennent de `residues/data/structure_components.yaml`, sourcé sur le dictionnaire des composants chimiques de la wwPDB. Un saut de numérotation y est signalé, jamais pris pour une rupture.
- **Extraction d'un peptide** (`io/peptide_extraction.py`, revue de Firas du 27 novembre 2026) : `analyze_chain` produit toujours un **rapport d'extraction** (`ExtractionReport`) ; `extract_peptides(structure, chaîne, journal=..., decisions=...)` en tire le ou les `Peptide`, ou lève `ExtractionRefusedError`, qui porte le rapport ; `structure_to_peptide` en est la forme stricte (un seul peptide, aucune décision). Chaque **jonction** C(i)–N(i+1) reçoit les résidus, les atomes utilisés (avec leur conformation), la distance mesurée, une catégorie (`continuous`, `numbering_discontinuity`, `break_declared`, `break_undeclared`, `continuity_unmeasurable`) et la décision appliquée. La liaison est présente si la distance ne dépasse pas le seuil de `residues/data/structure_geometry.yaml` (1,8 Å, valeur par défaut de `PPBuilder` de Biopython 1.88, validée par Firas ; longueur de référence 1,329 Å, Engh et Huber 1991), modifiable par `load_peptide_bond_criterion(valeur)`, la source devenant alors « valeur imposée ». Chaque **composant** hors squelette est classé : `terminal_cap`, `modified_residue` (atomes N, CA, C, absent de la bibliothèque), `side_chain_modification` (relié à la chaîne, directement ou par d'autres composants, par une liaison covalente **déclarée**, comme un glycane sur Asn), `unbound_ligand` (avec sa plus courte distance à la chaîne, sans seuil), `ion`, `water`. Les résidus déclarés non observés sont placés aux extrémités (`unobserved_termini`) ou dans une rupture ; un résidu déclaré qui n'a pas de place est une contradiction, refusée. Les atomes déclarés non observés (`REMARK 470`, `pdbx_unobs_or_zero_occ_atoms`) et les atomes de squelette absents sont rapportés à part. **Strict par défaut** : rupture, continuité non mesurable, extrémité non observée, résidu modifié, modification de chaîne latérale, liaison non représentable sont refusés ; les décisions explicites sont `split` (ruptures) et `observed_segment` (extrémités non observées), journalisées avec le rapport. Rien n'est jamais reconnecté ni reconstruit. Ligands, ions et eaux sont exclus du peptide et listés, jamais intégrés à la séquence.
- **Écart avec le plan, justifié** : la mesure C–N appartenait au contrôle géométrique de l'étape 4. L'étape 3 n'en avance que ce que l'extraction exige (une distance par jonction), parce qu'une séquence lue à travers une rupture serait fausse. Les autres contrôles de l'étape 4 (atomes manquants comparés au dictionnaire, stéréochimie, chevauchements) y restent.
- **Ne fait jamais** : corriger silencieusement une structure lue. Toute correction est un `StructureModifier` de `build`, journalisé.

### 5.8 `build` — construction et correction 3D

- **Interfaces** : `StructureBuilder.build(peptide, journal) -> Structure` ; `StructureModifier.apply(structure, journal) -> Structure`.
- **Erreurs** : `StructureError`, `ChemistryError`.
- **Garanties** : structure d'entrée jamais modifiée sur place ; extrémités conservées ; chaque transformation journalisée.
- **Ne fait jamais** : utiliser une nomenclature de champ de force.

### 5.9 `forcefields` — adaptateurs de champ de force

- **Interface** : `ForceFieldAdapter` avec `name`, `supported_engines`, `compatible_water_models`, `parametrize(peptide) -> ParameterSet`.
- **Erreurs** : `MissingParameterError` dès qu'un paramètre manque ; `ConfigurationError` pour un modèle d'eau incompatible.
- **Données** : `forcefields/data/` (correspondances de noms de résidus, fichiers à charger), chacune avec sa source.
- **Ne fait jamais** : écrire un fichier moteur, inventer un paramètre, substituer `NME` à `NHE`.

### 5.10 `engines` — fichiers d'entrée par moteur

- **Interface** : `EngineWriter` avec `engine`, `write_topology(structure, parameters, outdir) -> EngineFiles`, `solvate(files, spec, outdir, journal) -> EngineFiles`. Une méthode `write_protocol` (minimisation, équilibration, production) sera ajoutée à l'étape 9 du plan.
- **Erreurs** : `MissingParameterError` si `parameters.missing` n'est pas vide ; `ConfigurationError` si le couple champ de force / moteur n'est pas pris en charge ; `ExternalToolError`.
- **Garantie** : charge nette et nombre d'atomes relus dans les fichiers produits.
- **Ne fait jamais** : convertir implicitement depuis un autre moteur.

### 5.11 `pipeline` — orchestration

- **Interfaces** : `Stage[In, Out]` (`name`, `run(data, context) -> Out`), `StageContext`, `run_stage(stage, data, context, checks)`.
- **Rôle** : choisir les composants d'après `PrepConfig`, les enchaîner, appliquer les contrôles après chaque étape.
- **Ne fait jamais** : contenir de la chimie. Toute logique scientifique vit dans les modules des couches inférieures.

### 5.12 `cli` — interface en ligne de commande

- **Rôle** : traduire des arguments en appels à `pipeline` ou `io`, afficher les résultats, convertir `PepForgeError` en message lisible et code de sortie 1.
- **Commandes actuelles** : `pepforge parse --notation 1|3 SÉQUENCE` (analyse, affichage dans les deux codes) ; `pepforge convert --notation 1|3 SÉQUENCE` (réécriture dans l'autre code, seul le résultat sur la sortie standard, avertissements sur la sortie d'erreur) ; `pepforge config-template`. Sans `--notation`, `parse` et `convert` s'arrêtent sur un message qui dit quoi ajouter.

---

## 6. Flux du pipeline

### 6.1 Les trois modes d'entrée convergent

```mermaid
flowchart LR
    seq["Séquence<br/>Ac-RGDS-NH2"] -->|io.parse_sequence| pep[Peptide]
    pdb["PDB / mmCIF<br/>sortie de docking"] -->|io.read_structure| st0[Structure importée]
    smi["SMILES / MOL / SDF"] -->|io.read_fragment| frag[Fragment chimique]

    st0 -->|inventaire + correspondance| pep
    frag -->|forcefields : dérivation des paramètres| lib[(Bibliothèque locale<br/>de résidus)]
    lib -.-> pep

    pep -->|build.StructureBuilder| st1[Structure construite]
    st0 -->|build.StructureModifier<br/>nettoyage, atomes manquants, H| st1
```

Les trois modes aboutissent au même couple `Peptide` + `Structure`. Un fragment personnalisé n'entre pas directement dans la chaîne : il devient d'abord un résidu de la bibliothèque locale, avec ses paramètres et leur provenance, puis il est utilisé comme n'importe quel résidu.

### 6.2 Chaîne de préparation

```mermaid
flowchart TB
    cfg[/PrepConfig/] --> S1
    S1["1. Entrée<br/>io"] -->|Peptide ou Structure| C1{{contrôles 1}}
    C1 --> S2["2. Normalisation et protonation<br/>residues + build"]
    S2 -->|Peptide protoné| C2{{contrôles 2}}
    C2 --> S3["3. Construction ou correction 3D<br/>build"]
    S3 -->|Structure| C3{{contrôles 3}}
    C3 --> S4["4. Paramétrisation<br/>forcefields"]
    S4 -->|ParameterSet| C4{{contrôles 4}}
    C4 --> S5["5. Topologie dans le vide<br/>engines.write_topology"]
    S5 -->|EngineFiles| C5{{contrôles 5}}
    C5 --> S6["6. Solvatation, ions<br/>engines.solvate"]
    S6 -->|EngineFiles solvatés| C6{{contrôles 6}}
    C6 --> S7["7. Scripts de protocole<br/>engines.write_protocol"]
    S7 -->|EngineFiles + scripts| C7{{contrôles 7}}
    C7 --> S8["8. Rapport reproductible<br/>report"]
    C1 & C2 & C3 & C4 & C5 & C6 & C7 -. FAIL .-> STOP([QualityCheckError :<br/>arrêt + rapport partiel])
```

> **Écart avec l'ordre de la spécification (section 4 de [01-specification-fonctionnelle.md](01-specification-fonctionnelle.md))** : la construction 3D (étape 3 ici) précède la paramétrisation (étape 4). La structure est construite dans une nomenclature neutre, puis le champ de force est appliqué ; on peut ainsi changer de champ de force sans reconstruire. Le « choix du couple champ de force / moteur » de la spécification n'est pas une étape : il est fixé dans `PrepConfig` et vérifié avant le lancement (`require_complete()`).

> **Conversion entre moteurs** (ex. AMBER → GROMACS) : étape optionnelle 5 bis, distincte, avec ses propres contrôles (nombre d'atomes, charge nette, énergie en point unique comparés avant et après). Jamais implicite.

### 6.3 Exécution d'une étape

```mermaid
sequenceDiagram
    participant P as pipeline
    participant S as Stage
    participant J as Journal
    participant C as checks
    P->>S: run(data, context)
    S->>J: record(quoi, pourquoi, paramètres)
    S-->>P: sortie typée
    P->>C: run_checks(contrôles de l'étape, sortie)
    C-->>P: list[CheckResult]
    P->>J: record("stage_completed", avertissements)
    alt au moins un FAIL
        P-->>P: QualityCheckError (résultats conservés pour le rapport)
    else
        P-->>P: étape suivante
    end
```

### 6.4 Tableau des étapes

| # | Étape | Module | Entrée → sortie | Contrôles à la sortie (items de la checklist 04) |
|---|---|---|---|---|
| 1 | Entrée | `io` | texte ou fichier → `Peptide` / `Structure` | résidus reconnus, extrémités déclarées (sinon `WARN`), inventaire de la structure importée |
| 2 | Normalisation, protonation | `residues` (`build` n'intervient pas : aucune géométrie) | `Peptide` + pH → `ProtonationResult` (dont le `Peptide` protoné) | `PROTONATION_CHECKS` : états portés = états décidés, charge nette entière = somme des états, chaque pKa sourcé, état déduit = forme majoritaire ; `WARN` si \|pH − pKa\| < 1 ou si un état imposé est minoritaire |
| 3 | Construction / correction 3D | `build` | `Peptide` ou `Structure` → `Structure` | nombre de résidus, extrémités conservées, aucun atome manquant, hydrogènes présents, stéréochimie, conformations alternatives résolues |
| 4 | Paramétrisation | `forcefields` | `Peptide` → `ParameterSet` | `missing` vide, provenance renseignée, modèle d'eau compatible, fragments obtenus par analogie signalés |
| 5 | Topologie | `engines` | `Structure` + `ParameterSet` → `EngineFiles` | charge nette relue = charge théorique, nombre d'atomes, résidus terminaux (`NHE` et non `NME` pour `-CO-NH2` en AMBER), aucun avertissement de l'outil |
| 6 | Solvatation | `engines` | `EngineFiles` + `SolvationSpec` → `EngineFiles` | système neutre, concentration en sel, marge de boîte, aucun chevauchement |
| 7 | Protocole | `engines` | `EngineFiles` + protocole → scripts | équilibration NVT puis NPT, pas d'intégration cohérent avec les contraintes, électrostatique longue portée déclarée, répliques |
| 8 | Rapport | `report` | tout ce qui précède → rapport + YAML rejouable | versions d'outils, journal complet, objectif de simulation consigné |

---

## 7. Vérification : revue, contrôles et tests

Chaque étape est vérifiée à quatre niveaux. Les trois premiers sont automatiques et tournent dans la CI à chaque PR ; le quatrième combine automatisation et relecture scientifique.

```mermaid
flowchart LR
    U["1. Tests unitaires<br/>une fonction, un cas"] --> K["2. Tests de contrat<br/>une implémentation respecte son Protocol"]
    K --> I["3. Tests d'intégration<br/>plusieurs étapes enchaînées"]
    I --> V["4. Validation scientifique<br/>peptides de référence, outils réels"]
    A["Test d'architecture<br/>matrice de dépendances"] -.-> U
```

### 7.1 Niveaux de tests

| Niveau | Ce qui est vérifié | Où | Quand il tourne |
|---|---|---|---|
| **Architecture** | Imports conformes à la matrice, pas de cycle, tout module déclaré | `tests/test_architecture.py` | à chaque PR (CI) |
| **Unitaire** | Une fonction : cas nominal, cas limites, erreurs attendues | `tests/test_<module>.py` | à chaque PR (CI) |
| **Contrat** | Une implémentation respecte son `Protocol` et les garanties de la section 5 | `tests/contracts.py` (fonctions réutilisables), appelées par les tests de chaque implémentation | à chaque PR (CI) |
| **Intégration** | Enchaînement d'étapes via `run_stage`, avec composants factices ou réels | `tests/integration/` | à chaque PR ; les tests nécessitant AmberTools, GROMACS ou NAMD sont marqués `@pytest.mark.external` et sautés si l'outil est absent |
| **Validation scientifique** | Résultats comparés à des valeurs de référence établies indépendamment (à la main ou par un outil publié) | `tests/reference/` + relecture | avant chaque version, et à chaque ajout de résidu, de champ de force ou de moteur |

Fonctions de contrat disponibles dans [`tests/contracts.py`](../tests/contracts.py) :

- `assert_forcefield_contract(adapter, peptides)` : `ParameterSet` complet, `missing` vide, provenance présente ;
- `assert_engine_writer_contract(writer, structure, tmp_path)` : refus d'un `ParameterSet` incomplet ;
- `assert_builder_contract(builder, peptide)` : structure produite, extrémités conservées, construction journalisée ;
- `assert_stage_contract(stage, data, tmp_path)` : l'étape a un nom et s'exécute avec un contexte.

### 7.2 Peptides de référence

Jeu commun à toutes les étapes de validation, enrichi au fil du plan. Les valeurs attendues sont établies indépendamment de PepForge et documentées avec leur source dans `tests/reference/` avant d'être utilisées.

| Peptide | Ce qu'il exerce |
|---|---|
| `Ac-RGDS-NH2` | extrémités acétylée et amidée (`NHE` en AMBER), charges opposées Arg / Asp |
| `H-RGDS-OH` | extrémités libres, charges terminales |
| `Ac-RGDS-NHMe` | N-méthylamide : vérifie que `NME` n'est pas confondu avec `NHE` |
| Peptide avec His | choix du tautomère et de l'état de protonation selon le pH |
| Peptide avec Cys–Cys | pont disulfure : déclaré dans `core` et inventorié depuis l'étape 1 ; sa présence dans le `.psf` se vérifie à l'étape 7 |
| `Ac-RGDS-NH2` avec une position substituée en aza-β³ | **cas d'usage central** : greffage, paramètres dérivés, torsion N–N (section 10) |
| Le même peptide substitué en chaque position | mode lot, comparaison entre positions |
| Peptide portant deux résidus aza-β³ **adjacents** | charges de jonction entre deux résidus dérivés (section 12.1) |
| Dipeptide, et peptide d'une centaine de résidus | les deux extrémités de la gamme 2 à 100 résidus (section 12.5) |
| Peptide portant un chélateur, avec et sans métal | volet radiomarquage (sections 11.2 et 11.4) |
| Peptides de la thèse (ciblant DCC, TrkC, UNC5) | à fournir par Firas |

### 7.3 Contrôles d'exécution

Les contrôles de `checks` ne sont pas des tests : ils tournent à chaque préparation réelle, sur les données de l'utilisateur. Un contrôle couvre un item de [04-checklist-qualite.md](04-checklist-qualite.md) et l'indique dans `checklist_item`. Objectif : que chaque item automatisable de la checklist ait son contrôle. Les items non automatisables (relecture de la convergence de l'équilibration, par exemple) apparaissent dans le rapport comme « à vérifier manuellement ».

### 7.4 Revue de chaque PR

Le modèle [`.github/pull_request_template.md`](../.github/pull_request_template.md) impose, pour chaque PR :

1. l'étape du plan concernée et ses critères d'acceptation, cochés un à un ;
2. les modules touchés et la confirmation que la matrice de dépendances n'a pas changé (ou la justification de son changement, reportée dans ce document) ;
3. les tests ajoutés à chaque niveau ;
4. la provenance de toute nouvelle donnée scientifique ;
5. la CI verte : `ruff`, `mypy --strict`, `pytest` sur Python 3.11, 3.12 et 3.13.

---

## 8. Modifier un module isolément

| Je veux... | Je touche | Je ne touche pas | Je vérifie |
|---|---|---|---|
| Ajouter un acide aminé non standard | `residues/data/` (+ `forcefields/data/` pour ses paramètres) | le code de `io`, `build`, `engines` | test de charge nette et de paramètres manquants (règle 7) |
| Ajouter un champ de force | un nouvel adaptateur dans `forcefields/` | `engines` (s'il cible un moteur existant) | `assert_forcefield_contract` + validation sur les peptides de référence |
| Ajouter un moteur | un nouveau `EngineWriter` dans `engines/` | `forcefields`, `build` | `assert_engine_writer_contract` + comparaison avec un moteur existant |
| Changer d'outil de construction 3D | un nouveau `StructureBuilder` dans `build/` | `forcefields`, `engines` | `assert_builder_contract` + contrôles de l'étape 3 |
| Lire un nouveau format | `io/` | tout le reste | tests unitaires de lecture, aller-retour lecture / écriture |
| Ajouter un contrôle qualité | `checks/` et la déclaration de l'étape dans `pipeline` | les composants contrôlés | test du contrôle sur un cas `PASS` et un cas `FAIL` |
| Ajouter une commande | `cli/` | toute la logique | test de la commande (`tests/test_cli.py`) |

Si une modification oblige à toucher deux modules qui ne devraient pas se connaître, c'est le signe qu'il manque un type d'échange dans `core` ou un paramètre dans `PrepConfig` : ajouter ce type plutôt qu'un import de contournement.

---

## 9. Décisions ouvertes et compromis

### 9.1 Décisions scientifiques en attente (Firas)

> Mise à jour du 26 novembre 2026 : le descriptif de thèse (section 11) tranche le moteur, le système simulé et la variante pseudopeptidique. Les lignes barrées ci-dessous sont résolues.

| Décision | Bloque | Où elle s'enregistre |
|---|---|---|
| ~~Moteur principal~~ → **NAMD**, d'après le projet de thèse | étapes 7 et suivantes du plan | `PrepConfig.engine = namd` |
| ~~Champ de force principal~~ → **CHARMM36m + CGenFF**, confirmé par Firas le 26 novembre 2026 | étape 6 du plan | `PrepConfig.forcefield` |
| ~~Système simulé~~ → **les deux**, le complexe peptide–récepteur étant le cas principal | étapes 10 et 12 du plan | `PrepConfig.goal` |
| ~~Périmètre du radiomarquage~~ → **les deux** : peptide + chélateur sans métal, **et** peptide marqué (section 11.2) | volet chélateurs | `PrepConfig` (champ à ajouter) |
| Approche de paramétrisation du centre métallique (non lié, lié, atomes factices) et format retenu pour le système marqué (section 11.4) | étape 11F du plan | procédure de dérivation + `PrepConfig.forcefield` |
| Outil de construction 3D | étape 5 du plan | composant de `build` |
| Bibliothèque chimique pour SMILES / SDF (RDKit ou autre) | étape 11A du plan | dépendance optionnelle de `pyproject.toml` |
| Variante aza exacte : aza-β³ `-NH-N(R)-CH2-CO-` ou azapeptide α `-NH-N(R)-CO-` | tout le fil aza (section 10) | `ResidueTemplate` de la bibliothèque |
| Méthode de charges et niveau de théorie quantique | étape 11C du plan | procédure de dérivation, consignée dans le rapport |
| Logiciel de calcul quantique disponible — **non mentionné dans le projet de thèse, reste le principal inconnu** | étapes 11C et 11D du plan | `ExternalTool` dédié |
| ~~Séquences candidates~~ → **tout peptide de 2 à 100 résidus**, pas de liste figée | peptides de référence des tests | — |
| Moteurs : **NAMD d'abord, GROMACS ensuite** (même champ de force), **AMBER en dernier** (coût élevé, section 12.4) | étape 13 du plan | `PrepConfig.engine` |

### 9.2 Compromis d'architecture assumés

- **Construction 3D et paramétrisation séparées.** Des outils comme `tleap` font les deux à la fois. Ici, `build` produit une structure neutre et `engines` peut rappeler `tleap` pour écrire la topologie. Coût : une étape de plus. Bénéfice : changer de champ de force sans reconstruire, et tester chaque étape séparément.
- **`ParameterSet` décrit, il ne contient pas toutes les constantes.** Pour un champ de force publié, il référence les fichiers officiels à charger (`parameter_files`) plutôt que de recopier chaque constante dans des objets Python. Les constantes restent dans leurs fichiers d'origine, ce qui respecte la règle « pas d'invention de paramètres ».
- **`Protocol` plutôt que classes abstraites.** Une implémentation n'hérite de rien : elle est remplaçable par n'importe quel objet de même forme. Les tests de contrat compensent l'absence de vérification à l'héritage.
- **Contrôles séparés des composants.** Un même contrôle (charge nette, extrémités) s'applique à tous les moteurs et champs de force, au lieu d'être réécrit dans chaque implémentation.
- **`solvation` hors de `missing_decisions()`.** Les paramètres de solvatation ne sont exigés qu'au lancement de l'étape de solvatation, pour pouvoir utiliser PepForge sur les premières étapes sans les avoir fixés.

---

## 10. Cas d'usage central : substituer un acide aminé par un résidu aza-β³

Cette section a été ajoutée après la précision de Firas (26 novembre 2026) : la fonction centrale de PepForge est de **remplacer un acide aminé d'un peptide par le résidu pseudopeptidique aza-β³ correspondant, pour n'importe quel acide aminé et n'importe quelle position, puis de produire les fichiers de dynamique moléculaire**. Ce n'est donc pas un cas particulier du mode « structure chimique personnalisée » : c'est le chemin principal, et le reste de l'outil doit le servir.

### 10.1 Chaîne attendue

```mermaid
flowchart LR
    pep["Peptide de départ<br/>séquence ou structure 3D"] --> sub["Substitution<br/>position i → aza-β³"]
    lib[("Bibliothèque aza-β³<br/>un résidu par chaîne latérale<br/>× état de protonation")] --> sub
    sub --> st["Structure 3D du pseudopeptide"]
    st --> par["Paramétrisation<br/>résidu dérivé + champ de force du reste"]
    par --> files["Topologie, paramètres,<br/>coordonnées, scripts"]
    sub -.->|mode lot| scan["Scan de toutes les positions<br/>i = 1..n"]
```

### 10.2 Ce que la substitution change, du point de vue du modèle

> **À confirmer par Firas** : la variante visée est bien l'aza-β³, c'est-à-dire un squelette `-NH-N(R)-CH2-CO-`. L'azapeptide « α » classique (`-NH-N(R)-CO-`, semicarbazide, sans le CH2) est une autre famille. Les conséquences ci-dessous sont écrites pour l'aza-β³ ; elles changent si c'est l'autre variante.

| | Résidu α standard | Résidu aza-β³ |
|---|---|---|
| Squelette | `-NH-CαH(R)-CO-` | `-NH-N(R)-CH2-CO-` |
| Atomes de squelette | 3 (N, Cα, C) | 4 (N, N, CH2, C) |
| Liaison nouvelle | — | **N–N** (hydrazide) |
| Centre stéréogène | Cα (L ou D) | **disparaît** : l'azote substitué n'est pas un centre stéréogène stable |
| Chaîne latérale | portée par Cα | portée par l'azote |
| Torsions du squelette | φ, ψ, ω | φ, ψ, ω **plus** la torsion autour de N–N |

Trois conséquences pour l'architecture :

1. La substitution **allonge le squelette d'un atome**. Le modèle de données ne peut donc pas supposer qu'un résidu vaut trois atomes de squelette ; un contrôle doit vérifier la connectivité après greffage plutôt que de compter des atomes supposés.
2. Elle **supprime un centre stéréogène** : les contrôles de stéréochimie doivent attendre un Cα chiral sur les résidus α et aucun sur les résidus aza.
3. Elle **introduit une liaison N–N dans le squelette**, ce qui est le point dur ci-dessous.

### 10.3 Le point dur : aucun champ de force protéique ne couvre le squelette N–N

Les champs de force protéiques (ff14SB, ff19SB, CHARMM36m) décrivent le squelette des acides aminés α ; leurs torsions de squelette sont ajustées sur ce squelette. Un résidu aza-β³ n'y figure pas. Les champs de force généralistes (GAFF2, CGenFF) couvrent les molécules organiques, hydrazides compris, mais par analogie : leurs torsions autour d'une liaison N–N ne sont pas ajustées pour reproduire les préférences conformationnelles d'un pseudopeptide.

En conséquence, pour **chaque** résidu aza-β³ de la bibliothèque :

- les paramètres doivent être **dérivés** (voie GAFF2 via `antechamber` / `parmchk2` pour AMBER, CGenFF pour CHARMM) et non repris d'un champ de force publié ;
- les charges doivent être calculées sur un **fragment coiffé** (type `Ace-X-Nme`), puis ajustées pour que la charge nette du résidu soit entière et que les atomes de jonction restent compatibles avec les résidus voisins du champ de force protéique ;
- les **torsions du squelette, et en particulier celle autour de N–N, doivent être validées** contre des balayages dièdres relaxés en calcul quantique, avant toute simulation de production. Sans cette validation, les conformations obtenues ne sont pas interprétables.
- la **recherche bibliographique** (existe-t-il des paramètres publiés pour des aza-peptides ou des β-peptides réutilisables, et à quel niveau de théorie ?) est une tâche à part entière, préalable à l'étape de dérivation. Aucun jeu de paramètres n'est écrit dans le dépôt sans sa source ou sa procédure de dérivation complète.

### 10.4 Ce que l'architecture actuelle couvre déjà

| Besoin du cas d'usage | Déjà prévu |
|---|---|
| Représenter un peptide dont un résidu n'est pas standard | `core.Residue` a un indicateur `standard`, et `Peptide` n'impose aucune nomenclature |
| Remplacer un résidu dans une structure 3D | `build.StructureModifier` (mutation), avec structure immuable et journal obligatoire |
| Bibliothèque de résidus non standards versionnée, avec provenance | `residues/` et `residues/data/`, où chaque fichier déclare sa source |
| Refuser de simuler avec un paramètre manquant | `MissingParameterError`, `ParameterSet.missing` vide exigé, contrôle avant écriture moteur |
| Signaler un paramètre obtenu par analogie | `ParameterSource` et le niveau de confiance prévu à l'étape 11 du plan |
| Charge nette vérifiée | charge théorique (étape 2 du plan) comparée à la charge relue dans les fichiers produits (étape 7) |
| Entrée par SMILES / MOL / SDF pour un fragment | mode d'entrée 2.3 de la spécification, module `io` |

L'ossature est donc la bonne : **rien à reconstruire**. Ce qui manque est ci-dessous.

### 10.5 Ce qu'il faut ajouter

| # | À ajouter | Module | Pourquoi |
|---|---|---|---|
| A1 | Type `ResidueTemplate` : fragment chimique, atomes d'ancrage côté N et côté C, atomes à retirer lors du greffage, charge du fragment | `core` | C'est l'objet qu'échangent la bibliothèque, le greffage et la paramétrisation. Sans lui, `build` et `forcefields` devraient se connaître. |
| A2 | Bibliothèque aza-β³ : une entrée par chaîne latérale **et par état de protonation**, versionnée, avec la procédure de dérivation de chaque entrée | `residues/data/` | « Fonctionner sur tous les acides aminés » veut dire une vingtaine de chaînes latérales, dont sept ionisables, donc une trentaine d'entrées à dériver et valider. |
| A3 | Sous-pipeline de dérivation : fragment coiffé → calcul quantique → charges → `parmchk2` / CGenFF → fichier de paramètres + niveau de confiance | `forcefields` (nouveau sous-module `derivation`) | C'est l'étape la plus coûteuse et la plus sensible. Elle doit être reproductible, journalisée et rejouable à l'identique. |
| A4 | `MutateResidue` : greffage du template à la position i, reconstruction des liaisons de jonction, minimisation locale du voisinage, conservation du reste de la structure | `build` | Le cœur de la fonction demandée. |
| A5 | Mode « scan de positions » : même peptide, substitution en i pour i = 1..n, un système préparé par position | `pipeline` | C'est l'usage réel : comparer les positions entre elles, pas préparer un seul système. |
| A6 | Contrôles propres à la substitution : connectivité du squelette après greffage, longueur et angles de la liaison N–N dans une plage sourcée, absence de centre stéréogène attendu sur le résidu aza, charge nette entière du résidu, aucune torsion manquante autour de N–N, chaîne latérale intacte | `checks` | Les contrôles existants supposent un squelette α. |
| A7 | Validation des torsions contre des balayages dièdres quantiques, avec les courbes conservées dans `tests/reference/` | tests + `docs` | Sans elle, les paramètres dérivés ne sont pas défendables dans une thèse. |
| A8 | Journal et rapport étendus : pour chaque résidu dérivé, méthode de charges, niveau de théorie, pénalités d'analogie, version des outils | `report` | Traçabilité exigée par la règle 3 du projet. |

Ces ajouts ne changent **aucune** règle de dépendance de la section 3 : `ResidueTemplate` vit dans `core`, la dérivation dans `forcefields`, le greffage dans `build`, le lot dans `pipeline`. C'est la vérification que l'architecture tient pour ce cas d'usage.

### 10.6 Décisions à prendre avant de coder ce fil

| Décision | Pourquoi elle bloque |
|---|---|
| Variante exacte : aza-β³ (`-NH-N(R)-CH2-CO-`) ou azapeptide α (`-NH-N(R)-CO-`) | Change le template, les contrôles et les paramètres à dériver |
| Voie de paramétrisation : AMBER + GAFF2 ou CHARMM + CGenFF | Change tout le sous-module de dérivation |
| Méthode de charges et niveau de théorie (AM1-BCC, RESP sur HF/6-31G\*, autre) | Doit être fixée une fois pour toute la bibliothèque, sinon les résidus ne sont pas comparables entre eux |
| Logiciel de calcul quantique disponible | Conditionne la faisabilité de A3 et A7 |
| Substitution d'un seul résidu à la fois, ou plusieurs simultanément | Change `MutateResidue` et le mode lot |

---

## 11. Contraintes venues du projet de thèse

Section ajoutée le 26 novembre 2026 après lecture du descriptif de thèse de Firas (« Conception, synthèse, activité et validation in vitro de pseudopeptides théranostiques contenant des acides aminés aza-β³, spécifiques aux récepteurs à dépendance »). Les éléments ci-dessous viennent de ce document et remplacent plusieurs suppositions laissées ouvertes.

### 11.1 Ce que le document tranche

| Point | Ce que dit le projet de thèse | Conséquence pour PepForge |
|---|---|---|
| Variante pseudopeptidique | acides aminés **aza-β³** | confirmée par Firas ; la section 10 s'applique telle quelle |
| Champ de force | déduit du moteur NAMD | **CHARMM36m + CGenFF**, confirmé par Firas le 26 novembre 2026 ; fichiers détaillés en 11.3 |
| Moteur de dynamique moléculaire | **NAMD**, cité pour l'affinement des modèles, la dynamique moléculaire et la minimisation d'énergie, avec **VMD** | le premier `EngineWriter` est NAMD (`psfgen` : `.psf`, `.pdb`, plus `.prm`/`.rtf`/`.str`) |
| Visualisation | PyMOL, VMD | hors périmètre de PepForge, mais les fichiers produits doivent s'ouvrir dans ces outils |
| Système simulé | docking (AutoDock4, Vina, HADDOCK, PatchDock, FireDock) **puis** dynamique moléculaire pour évaluer affinité et spécificité | les **deux** scénarios sont nécessaires, et le **complexe peptide–récepteur** est le cas principal, pas une extension |
| Peptides de départ | régions de la **nétrine-1** et de la **NT-3**, séquences et structures tirées de la PDB et d'UniProt | les peptides de référence doivent venir de ces séquences ; à demander à Firas (deux projets de master déjà réalisés) |
| Récepteurs | DCC, TrkC, UNC5, récepteur de la nétrine-1 | structures à récupérer dans la PDB, à préparer avec les mêmes contrôles que le peptide |
| Structures expérimentales | RMN (CARA, ARIA), validation SAVES 6.0 | `io` doit lire un **ensemble RMN** (PDB multi-modèles) et pas seulement un modèle unique |

### 11.2 Ce que le document ajoute au périmètre : chélateurs et radiométaux

Les peptides de la thèse ne sont pas de simples peptides : ils portent un **chélateur** greffé à l'extrémité N-terminale ou C-terminale amide, et sont **radiomarqués**.

- Chélateurs prévus : DOTA (sous forme protégée DOTA(tBu)₃), N4 (6-carboxy-1,4,8,11-tétraazaundécane, (Boc)₄-N4), HYNIC (acide 6-hydrazinonicotinique, Boc-HYNIC).
- Radiométaux prévus : **gallium, rhénium, lutécium, technétium-99m**.
- Des groupes fonctionnels (DOTA, maléimide ou autres) peuvent être greffés sur la chaîne principale **ou sur une chaîne latérale**.

Conséquences, à trancher avec Firas avant d'écrire du code pour ce volet :

1. Un chélateur est, du point de vue de l'architecture, un **`ResidueTemplate` terminal ou latéral** : même mécanique de greffage et de dérivation que pour un résidu aza-β³ (section 10.5). L'architecture le couvre sans changement.
2. Le **complexe métallique est un problème distinct et plus difficile**. Les champs de force classiques ne décrivent pas les métaux de transition de façon transférable ; paramétrer un complexe Tc, Re, Lu ou Ga demande une approche dédiée (modèle de liaison, modèle de points cationiques factices, ou contraintes), chacune avec ses limites. Rien ne sera écrit sans source ni procédure.
3. **Périmètre tranché par Firas le 26 novembre 2026 : les deux cas.** PepForge doit préparer le peptide **portant son chélateur sans métal** et le peptide **marqué, métal compris**. Le premier cas est une extension directe de la section 10 ; le second est traité en 11.4, parce qu'il pose un problème de nature différente. C'est cohérent avec le projet de thèse, qui teste les peptides « marqués ou non ».

### 11.3 Les fichiers CHARMM et NAMD à produire

Le couple retenu est **CHARMM36m + CGenFF, moteur NAMD** (confirmé par Firas le 26 novembre 2026). Voici ce que PepForge doit produire, et qui le produit. La règle 2 du projet s'applique : l'adaptateur de champ de force décrit les paramètres, le module `engines` écrit les fichiers.

| Fichier | Ce qu'il contient | Produit par |
|---|---|---|
| `.rtf` (*residue topology file*) | pour chaque résidu : noms et types d'atomes, charges partielles, connectivité, et les *patches* (`PRES`) de liaison et d'extrémité | `forcefields` décrit, `engines` écrit |
| `.prm` ou `.par` (*parameter file*) | paramètres liés et non liés **par type d'atome** : liaisons, angles, dièdres, impropres, termes de Lennard-Jones | fourni par CHARMM36m pour les résidus standards ; complété pour les résidus dérivés |
| `.str` (*stream file*) | un fichier qui porte **à la fois** une section de topologie et une section de paramètres, pour une molécule ou un résidu donné. C'est la sortie de CGenFF pour un fragment personnalisé | `forcefields.derivation` |
| `.psf` (*protein structure file*) | la description du système assemblé : tous les atomes, leurs types, leurs charges, et les listes de liaisons, angles, dièdres et impropres | `engines` via `psfgen` |
| `.pdb` ou `.crd` | les coordonnées | `engines` |
| fichier de configuration NAMD | le script de la simulation : minimisation, équilibration, production, traitement de l'électrostatique, rayons de coupure | `engines.write_protocol` |

Deux pièges à éviter, inscrits comme contrôles :

1. **Un `.str` n'est pas un fichier de coordonnées.** La spécification fonctionnelle le signalait déjà ; le contrôle refuse un `.str` présenté comme tel.
2. **La sortie de CGenFF décrit une molécule isolée, pas un résidu de chaîne.** Pour intégrer un résidu aza-β³ dans un peptide, il faut en tirer une entrée `RESI` avec ses atomes de liaison, ou un `PRES`, puis **ajuster les charges pour que le résidu porte une charge entière et que ses atomes de jonction restent compatibles avec les résidus CHARMM36m voisins**. C'est le cœur de l'étape 11C, et c'est là que se joue la qualité du résultat. Les **pénalités CGenFF** (pénalité de charge et pénalité de paramètre) sont conservées et reportées dans le rapport : au-delà d'un seuil à fixer avec Firas, elles imposent une dérivation quantique.

### 11.4 Le centre métallique : un problème de nature différente

Le périmètre inclut le peptide marqué (gallium, rhénium, lutécium, technétium-99m). Il faut dire clairement pourquoi c'est plus difficile que le reste : les champs de force classiques décrivent mal les métaux de transition, parce que la coordination est directionnelle et partiellement covalente, alors que ces champs de force reposent sur des charges ponctuelles et des termes non liés isotropes.

Trois familles d'approches existent, et **le choix doit être documenté avec sa source** ; aucune n'est retenue par défaut :

| Approche | Principe | Limite principale |
|---|---|---|
| Modèle non lié | le métal est un ion ponctuel, charge et termes de Lennard-Jones | ne tient pas la géométrie de coordination : le métal peut quitter le chélateur |
| Modèle lié | liaisons explicites métal–ligand, paramètres ajustés sur un calcul quantique du complexe | demande un calcul quantique par complexe, et interdit tout échange de ligand |
| Atomes factices cationiques | des pseudo-atomes autour du métal imposent la géométrie | paramétrisation délicate, et dépend de la géométrie de coordination supposée |

Deux points d'honnêteté à porter dans la thèse :

- Le **technétium et le rhénium** sous forme de cœurs oxo (typiquement rencontrés avec HYNIC et des co-ligands) sont les cas les plus difficiles : la nature exacte du complexe formé lors du radiomarquage n'est pas toujours établie. Si elle ne l'est pas, le système simulé repose sur une hypothèse structurale qu'il faut énoncer comme telle. Le gallium et le lutécium, cations durs dans un DOTA, sont plus traitables.
- **Tension de format à surveiller.** L'outillage de construction des modèles liés métal–ligand est historiquement plus développé du côté AMBER que du côté CHARMM. Or NAMD lit aussi les topologies AMBER. Si un modèle lié s'avère nécessaire, il y a deux voies acceptables, et une interdite :
  - dériver ou trouver des paramètres **compatibles CHARMM**, et rester en `.psf` / `.prm` ;
  - préparer **tout le système marqué** au format AMBER, lu par NAMD, en traitant cela comme une préparation distincte et vérifiée ;
  - **interdit** : mélanger un champ de force protéique CHARMM et des paramètres métalliques AMBER dans un même système. C'est la règle 4 du projet, et un contrôle doit le refuser.

### 11.5 Cohérence avec le reste de l'architecture

Aucune règle de dépendance de la section 3 ne change :

- le chélateur et le résidu aza-β³ passent tous deux par `core.ResidueTemplate`, la bibliothèque de `residues/` et la dérivation de `forcefields` ;
- NAMD devient une implémentation de `EngineWriter` dans `engines`, sans effet sur `forcefields` ni sur `build` ;
- l'ensemble RMN multi-modèles est une capacité de lecture de `io`, avec un contrôle qui impose le choix explicite du modèle à préparer ;
- le complexe peptide–récepteur reste l'étape 12 du plan, mais devient prioritaire au même titre que le peptide libre.

---

## 12. Portée confirmée le 26 novembre 2026

Firas a précisé quatre points qui élargissent la portée. Chacun est traité ci-dessous avec son coût réel.

### 12.1 Substitutions multiples, en nombre et en types quelconques

Une préparation doit accepter **une, deux ou plus** substitutions aza-β³, sur un **même type** d'acide aminé ou sur **plusieurs types différents**, dans un peptide de **2 à 100 résidus**.

Conséquences :

- La substitution est une opération **répétable**, pas un cas unique : `MutateResidue` s'applique n fois, chaque application étant journalisée séparément.
- Le cas **deux résidus aza-β³ adjacents** est particulier : les deux atomes de jonction qui se font face appartiennent tous deux à un résidu dérivé, et non à un résidu CHARMM36m. Il demande son propre jeu de charges de jonction et son propre test.
- Le mode lot de l'étape 11E ne se limite pas à « une position à la fois » : il doit pouvoir énumérer des **combinaisons** de positions. Le nombre de combinaisons explose avec la taille du peptide, donc l'utilisateur déclare explicitement les combinaisons voulues ; PepForge ne les énumère jamais toutes de lui-même.

### 12.2 Deux voies de paramétrisation, et pourquoi la bibliothèque reste la voie principale

Firas souhaite faire tourner CGenFF sur **la molécule entière** modifiée. C'est faisable, mais ce n'est pas interchangeable avec la voie par bibliothèque de résidus. Voici la comparaison honnête.

| | Voie A — bibliothèque de résidus | Voie B — molécule entière dans CGenFF |
|---|---|---|
| Ce qui passe dans CGenFF | le fragment aza-β³ coiffé, une fois par chaîne latérale | le peptide modifié complet, à chaque nouveau peptide |
| Résidus standards | gardent **CHARMM36m**, avec ses torsions de squelette et sa correction **CMAP** ajustées pour les peptides | reçoivent des paramètres **génériques** CGenFF, sans CMAP |
| Réutilisation | une dérivation sert à tous les peptides et toutes les positions | aucune : chaque peptide, chaque position, chaque combinaison repart de zéro |
| Taille praticable | toute la gamme 2 à 100 résidus | petites molécules seulement ; la limite en nombre d'atomes de la version utilisée est à vérifier, et 100 résidus font plus d'un millier d'atomes |
| Jonctions | travail de charges de jonction à faire, une fois par chaîne latérale | pas de jonction à traiter : CGenFF voit la connectivité réelle |
| Pénalités | calculées sur le fragment, hors contexte de chaîne | calculées dans le contexte réel |

**Recommandation : la voie A est la voie principale, la voie B sert de contrôle croisé.** Deux raisons décisives :

1. La voie B **abandonne CHARMM36m pour les résidus standards**. Pour un peptide de 100 résidus, cela revient à décrire un petit peptide structuré avec un champ de force généraliste, sans la correction CMAP. C'est un recul, pas un progrès.
2. La voie B ne tient pas à l'échelle demandée : chaque combinaison de substitutions exigerait une paramétrisation complète.

En revanche, la voie B a une vraie valeur comme **vérification** : sur un peptide court (deux à cinq résidus) portant une substitution, on compare les charges et les paramètres obtenus par la bibliothèque à ceux obtenus sur la molécule entière. Un écart important révèle une erreur dans le traitement des jonctions. C'est un test de validation à part entière, inscrit à l'étape 11C.

**Limite à écrire dans la thèse** : un résidu aza-β³ ne bénéficie d'aucune correction CMAP, puisque celle-ci est ajustée pour le squelette des acides aminés α. C'est une raison de plus pour valider les torsions du squelette, y compris celle autour de N–N, contre des balayages quantiques (étape 11D).

### 12.3 « Un métal, quelle que soit sa nature »

Le **mécanisme** est générique, et c'est ce que l'architecture garantit : un centre métallique est décrit par un enregistrement de données (élément, état d'oxydation, charge, géométrie de coordination attendue, atomes donneurs du chélateur, approche de paramétrisation retenue), au même titre qu'un `ResidueTemplate`. Ajouter un métal, c'est ajouter une entrée de données, pas du code.

En revanche, et il faut être clair : **les paramètres, eux, ne sont pas génériques**. Pour chaque métal et chaque chélateur, il faut soit des paramètres publiés et cités, soit une dérivation documentée. Un métal sans entrée validée lève `MissingParameterError` ; PepForge ne fabrique jamais un jeu de paramètres métalliques par analogie silencieuse. Les limites propres à chaque approche restent celles de la section 11.4.

### 12.4 Trois moteurs : ce qui est peu coûteux et ce qui l'est beaucoup

Firas veut travailler sur NAMD avec CHARMM, mais pouvoir aussi produire des fichiers AMBER et GROMACS, et choisir selon la méthode. L'architecture est faite pour ça, mais les deux situations n'ont pas le même coût, et il ne faut pas les confondre.

| Situation | Ce que ça demande | Coût |
|---|---|---|
| **Même champ de force, autre moteur** — par exemple CHARMM36m exécuté sous GROMACS | un `EngineWriter` GROMACS, ou une conversion de format traitée comme une étape explicite et vérifiée | modéré : un générateur par moteur, les paramètres ne changent pas |
| **Autre champ de force** — par exemple AMBER ff19SB avec GAFF2 | un nouvel `EngineWriter`, un nouvel adaptateur, **et la dérivation complète de toute la bibliothèque aza-β³ une seconde fois** (GAFF2 via `antechamber` et `parmchk2`), plus la revalidation des torsions | élevé : c'est le travail des étapes 11C à 11E refait en entier |

Ordre de priorité qui en découle, inscrit au plan :

1. **CHARMM36m + NAMD** : la voie de la thèse, développée d'abord.
2. **GROMACS avec le même champ de force** : ensuite, comme conversion vérifiée (nombre d'atomes, charge nette et énergie en point unique comparés avant et après).
3. **AMBER ff19SB + GAFF2** : en dernier, parce que cela double la dérivation et la validation.

Interdit, et contrôlé : mélanger des paramètres de champs de force différents dans un même système, quelle que soit la combinaison de moteurs (règle 4 du projet).

### 12.5 Taille des peptides : de 2 à 100 résidus

Il n'y a pas de liste figée de séquences candidates : tout peptide est admis dans cette gamme. Conséquences :

- Les **peptides de référence** des tests restent des cas illustratifs choisis pour ce qu'ils exercent (section 7.2), et non les séquences de la thèse.
- Les tests doivent couvrir les deux extrémités de la gamme : un dipeptide et un peptide d'une centaine de résidus, ce dernier avec plusieurs substitutions.
- À 100 résidus, on prépare un petit domaine structuré : la correction CMAP de CHARMM36m et la qualité des torsions de squelette comptent davantage, ce qui renforce la conclusion de la section 12.2.

---

## 13. Intégrer un résidu aza-β³ dans CHARMM36m

Firas a posé le point juste, et c'est le vrai risque technique du projet : **paramétrer le fragment aza-β³ séparément ne garantit pas que son insertion dans CHARMM36m soit correcte**. Un aza-β³ n'est pas un acide aminé standard à chaîne latérale modifiée : il modifie le **squelette** lui-même. Tout ce qui, dans CHARMM36m, est ajusté pour le squelette des acides aminés α devient donc une question ouverte à cet endroit de la chaîne.

Cette section décrit la stratégie de jonction, les termes à vérifier un par un, le traitement de CMAP, et le protocole d'acceptation. Elle conditionne les étapes 11C et 11D du plan.

### 13.1 Les cinq familles de termes à vérifier

| Famille | Ce qui est en cause | Comment on le vérifie |
|---|---|---|
| **Liaisons de jonction** | les liaisons entre le résidu aza et ses voisins : `C(i-1)–N(i)` et `C(i)–N(i+1)` | types d'atomes des deux côtés présents dans `.prm` ; longueur d'équilibre comparée à l'optimisation quantique du tripeptide |
| **Angles du squelette** | tout angle dont le sommet ou une branche appartient au résidu aza, y compris ceux qui traversent la jonction | énumération exhaustive depuis le `.psf` produit, puis recherche de chaque triplet de types dans `.prm` |
| **Dièdres du squelette** | les dièdres qui traversent la jonction, l'amide `ω`, et les torsions **nouvelles** propres à l'aza-β³ : autour de `N–N`, de `N–CH2` et de `CH2–C(=O)` | balayages dièdres relaxés quantiques sur un tripeptide `Ac-X-aza-Y-NMe`, comparés au champ de force |
| **Impropres** | planéité du carbonyle, et surtout le caractère **pyramidal ou plan de l'azote substitué** : c'est lui qui gouverne la barrière d'inversion, donc la disparition du centre stéréogène | barrière d'inversion quantique comparée à celle du champ de force |
| **Non liées** | types de Lennard-Jones issus de CGenFF face à ceux de CHARMM36m, les **paramètres 1-4 explicites** de CHARMM (`eps14`, `rmin14`) et les éventuels termes `NBFIX` | les deux jeux partagent la même forme fonctionnelle et ont été développés pour être combinables, ce qui rend la voie viable ; les paramètres 1-4 du résidu dérivé doivent être déclarés explicitement là où CHARMM36m le fait |

Aucune de ces vérifications n'est facultative, et chacune devient un contrôle de `checks` plutôt qu'une relecture humaine : l'énumération se fait sur le `.psf` réellement produit, pas sur ce qu'on croit avoir écrit.

### 13.2 Stratégie de jonction : des patches explicites

Dans CHARMM, le lien entre résidus consécutifs est porté par des *patches* (`PRES`). Un résidu aza-β³ ne peut pas réutiliser le patch peptidique standard, puisque les atomes de jonction et leurs types diffèrent. Pour **chaque** type de résidu aza-β³, la bibliothèque doit donc fournir cinq patches :

| Patch | Jonction couverte |
|---|---|
| `standard → aza` | résidu α en i-1, résidu aza en i |
| `aza → standard` | résidu aza en i, résidu α en i+1 |
| `aza → aza` | deux résidus aza adjacents (cas signalé en 12.1) |
| `N-ter aza` | résidu aza en première position, avec l'extrémité choisie |
| `C-ter aza` | résidu aza en dernière position, avec l'extrémité choisie |

Deux règles qui en découlent, et qui deviennent des contrôles :

1. **Neutralité de charge à la jonction.** Après application du patch, la somme des charges partielles sur la jonction doit rester telle que le système entier porte la charge entière attendue. Un écart, même de 0,01 e, signale une erreur de redistribution.
2. **Aucune jonction sans patch.** Une combinaison de voisins pour laquelle aucun patch n'existe lève `MissingParameterError`. Elle n'est jamais rattrapée par le patch peptidique standard.

### 13.3 CMAP : ce qui s'applique, ce qui ne s'applique pas

C'est le point le plus délicat, et il faut être précis, parce que CMAP ne se comporte pas comme les autres termes.

Dans CHARMM36m, CMAP est une correction bidimensionnelle sur le couple `φ`/`ψ` du **squelette protéique**, déclarée résidu par résidu dans la topologie, et dont la définition **référence des atomes des résidus voisins** (le carbonyle du précédent, l'azote du suivant). D'où trois constats distincts :

1. **Le résidu aza-β³ ne reçoit aucun CMAP.** Sa suite d'atomes de squelette n'est pas celle d'un acide aminé α : les dièdres `φ`/`ψ` au sens de CMAP n'y sont pas définis, et il porte en plus une torsion que le squelette α ne connaît pas. Aucune grille existante ne lui est applicable. Ce n'est pas un réglage à trouver, c'est une absence.
2. **Les résidus voisins gardent leur CMAP, mais celui-ci traverse la jonction.** Le CMAP du résidu i+1 référence le carbonyle du résidu i, c'est-à-dire un atome du résidu aza. Deux choses à vérifier séparément : que la référence **se résout** vers un atome existant et typé (sinon `psfgen` échoue ou, pire, omet le terme silencieusement), et que la valeur de la correction reste **pertinente**, ce qui ne se décide que par comparaison quantique.
3. **La question « une correction complémentaire est-elle nécessaire ? » est ouverte, et son coût est très asymétrique.**

Trois options, avec leur coût réel :

| Option | Ce qu'elle demande | Quand la retenir |
|---|---|---|
| **A. Pas de CMAP sur le résidu aza, voisins inchangés, limite documentée** | les balayages dièdres de l'étape 11D, pour mesurer et consigner l'écart | voie par défaut : c'est la seule réaliste à l'échelle d'une bibliothèque d'une trentaine d'entrées |
| **B. Torsions du résidu aza réajustées pour compenser** | un ajustement des termes dièdres sur les profils quantiques, par type de résidu | si l'écart mesuré en A dépasse le seuil convenu |
| **C. Dérivation d'une grille de correction propre à l'aza-β³** | un balayage bidimensionnel quantique par type de résidu, soit des centaines de points chacun, puis l'implémentation du terme | hors de portée pour l'ensemble de la bibliothèque ; éventuellement sur **un** résidu de référence, comme contribution méthodologique à part entière |

**Recommandation : commencer par A, mesurer, puis passer à B pour les résidus dont l'écart le justifie.** L'option C ne se décide qu'au vu des chiffres de A, et vaut mieux comme résultat de thèse assumé que comme prérequis à l'outil. Dans tous les cas, **l'absence de CMAP sur le résidu modifié est inscrite dans le rapport de chaque préparation concernée**, pas seulement dans la documentation.

### 13.4 Protocole d'acceptation

Pour chaque type de résidu aza-β³ entrant dans la bibliothèque :

1. Construire le tripeptide de référence `Ac-X-aza-Y-NMe`, avec `X` et `Y` fixés une fois pour toute la bibliothèque afin que les entrées soient comparables.
2. Optimiser sa géométrie en quantique, au niveau de théorie unique retenu pour tout le projet.
3. Énumérer, depuis le `.psf` produit par `psfgen`, **tous** les termes liés et impropres impliquant le résidu aza et ses deux jonctions, puis vérifier que chacun existe dans les fichiers de paramètres. Un terme manquant lève `MissingParameterError`.
4. Lancer les balayages dièdres relaxés sur les torsions de la troisième ligne du tableau 13.1, et comparer les profils champ de force et quantique.
5. Comparer les géométries d'équilibre, la barrière d'inversion de l'azote substitué et la charge nette.
6. Consigner tous les écarts, les seuils appliqués et leur justification dans `tests/reference/`, puis dans le rapport.

Un résidu qui ne passe pas ce protocole n'entre pas dans la bibliothèque : il reste marqué comme non validé et son usage lève une erreur.

### 13.5 Effet sur l'architecture

Rien de tout cela ne change les règles de dépendance de la section 3 :

- les patches de jonction sont des **données** de `residues/data/` et de `forcefields/data/`, avec leur provenance ;
- l'énumération des termes depuis le `.psf` est un **contrôle** de `checks`, appliqué à la sortie de l'étape 5 du pipeline ;
- les balayages quantiques passent par `pepforge.tools`, comme tout exécutable externe ;
- la mention de l'absence de CMAP est une entrée de `report`.

---

## 14. Métaux : données livrées et ajout manuel

Firas demande que l'outil porte les paramètres de chaque métal **et** permette d'en ajouter manuellement. Les deux, et c'est le bon compromis.

### 14.1 Un métal est une entrée de données

Un centre métallique est décrit dans `forcefields/data/metals/<élément>-<chélateur>.yaml` :

| Champ | Contenu |
|---|---|
| `element`, `oxidation_state`, `charge` | identité du centre |
| `chelator` | chélateur auquel l'entrée s'applique (DOTA, N4, HYNIC...) |
| `coordination` | géométrie attendue et atomes donneurs |
| `approach` | `nonbonded`, `bonded` ou `cationic_dummy` (section 11.4) |
| `parameters` | les paramètres eux-mêmes, ou le chemin du fichier à charger |
| `source` | **obligatoire** : publication citée, ou procédure de dérivation complète |
| `confidence` | niveau de confiance et limites connues, repris dans le rapport |

Ajouter un métal ne demande donc **aucun code** : une entrée de données et son test. C'est ce qui rend la promesse « quel que soit le métal » tenable côté mécanisme.

### 14.2 Ce qui est livré, et ce qui ne peut pas l'être

PepForge livrera les entrées pour lesquelles il existe des paramètres publiés et citables, en commençant par les couples de la thèse : gallium et lutécium dans un DOTA, puis technétium et rhénium, plus difficiles (section 11.4). Chaque entrée livrée porte sa source.

Ce qui ne peut pas être livré : un jeu de paramètres pour un couple métal–chélateur qui n'a pas été publié ni dérivé. Dans ce cas, `MissingParameterError`. **L'outil ne complète jamais un métal par analogie silencieuse** : c'est la règle 2 du projet, et c'est aussi ce qui protège la thèse.

### 14.3 L'ajout manuel, avec sa traçabilité

Un utilisateur peut fournir sa propre entrée, par un fichier déclaré dans la configuration. Deux garde-fous, non négociables :

1. **`source` est obligatoire.** Une entrée sans source est refusée, pas acceptée avec un avertissement. Sans cela, la traçabilité du rapport ne vaut plus rien.
2. **Les mêmes contrôles s'appliquent** qu'à une entrée livrée : charge totale du complexe conforme, géométrie de coordination vérifiée après construction, aucun paramètre manquant, et interdiction de mélanger des champs de force (section 11.4).

L'entrée manuelle est recopiée intégralement dans le rapport, avec sa source, pour qu'une préparation reste rejouable même si le fichier d'origine a changé entre-temps.

---

## 15. Calcul quantique : recommandation et interface

Firas demande quel logiciel retenir. Voici la recommandation, puis la façon dont l'outil s'en abstrait.

### 15.1 Ce dont le projet a besoin

| Usage | Exigence |
|---|---|
| Optimisation de géométrie des fragments et tripeptides de référence | courant, aucun besoin exotique |
| Balayages dièdres relaxés, y compris autour de `N–N` (étapes 11D et 13.4) | nombreux calculs, donc **le coût et l'automatisation comptent plus que la fonctionnalité** |
| Barrière d'inversion de l'azote substitué | optimisation sous contrainte |
| Charges, si l'on va au-delà de CGenFF | dépend du protocole retenu, qui doit être fixé une fois pour toute la bibliothèque |

### 15.2 Recommandation

**ORCA en premier, Psi4 en second, Gaussian seulement si le laboratoire en a déjà la licence.**

| Logiciel | Pour | Contre |
|---|---|---|
| **ORCA** | gratuit pour l'usage académique sur enregistrement, très actif, bon rapport performance/précision en DFT et MP2, parallélisation correcte | pas libre au sens des sources ; conditions d'usage à relire |
| **Psi4** | logiciel libre, **interface Python** qui s'automatise parfaitement depuis PepForge, idéale pour les tests reproductibles | moins rapide qu'ORCA sur les gros calculs |
| **Gaussian** | référence historique, beaucoup de protocoles publiés le citent, utile si l'on veut reproduire exactement une procédure publiée | licence commerciale coûteuse ; aucun intérêt à l'acheter pour ce projet s'il n'est pas déjà disponible |

**Prendre en charge ORCA et Psi4 derrière la même interface**, ce qui coûte peu et évite d'être prisonnier d'un choix : si le laboratoire obtient Gaussian plus tard, c'est une troisième implémentation, pas une réécriture.

Le point qui compte davantage que le logiciel : **le niveau de théorie et le protocole sont fixés une fois pour toute la bibliothèque**, et consignés dans le rapport. Deux entrées dérivées à des niveaux différents ne sont pas comparables entre elles, et une bibliothèque hétérogène serait critiquable en soutenance. Le choix du niveau se fait au vu de la bibliographie de l'étape préalable au fil aza, pas avant.

### 15.3 Interface

Un protocole `QuantumTool` dans `pepforge.tools`, avec une implémentation par logiciel :

- `optimize(fragment, level) -> Geometry`
- `scan_dihedral(fragment, atoms, range, level) -> Profile`
- `charges(fragment, method, level) -> Charges`

Chaque appel passe par `ExternalTool` : version du logiciel enregistrée, sortie analysée pour détecter un calcul non convergé plutôt que de se fier au code de retour. Un calcul non convergé lève `ExternalToolError` et n'entre jamais dans une dérivation.

### 15.4 Ajouter concrètement ORCA et Psi4 au projet

Quatre choses à faire, dans cet ordre, et une seule demande du code.

**1. Installer les logiciels.** Ils ne s'installent pas par `pip install -e ".[dev]"`, parce que ce ne sont pas des dépendances Python de PepForge mais des exécutables externes, exactement comme `psfgen` ou `gmx`.

- **ORCA** : compte académique gratuit à créer sur le portail des développeurs, puis téléchargement de l'archive binaire pour la plateforme, extraction, et le dossier `bin` ajouté au `PATH`. Rien à compiler. La licence académique interdit la redistribution : l'archive n'entre jamais dans le dépôt.
- **Psi4** : installation par `conda` depuis son canal, dans un environnement dédié. C'est le plus simple des deux, et le seul qui s'importe aussi comme module Python, ce qui le rend idéal pour les tests reproductibles.

Aucun des deux ne devient une dépendance déclarée du paquet : PepForge les détecte sur le `PATH` au moment où il en a besoin, et lève `ExternalToolError` s'ils manquent.

**Installer ORCA sous Windows, pas à pas.** Compte à créer sur le forum ORCA (voie académique) ou dans l'espace de téléchargement de FACCTs, licence à accepter, puis téléchargement d'une archive `.zip` qui **contient un installateur** `.exe` : ce n'est pas une archive de binaires à poser à la main. L'archive est à décompresser, l'exécutable à lancer, et il faut choisir l'installation **« Custom » ou « Full »**, parce que l'installation « Typical » ne copie que les modules **série** et laisse de côté la partie parallèle. L'installateur écrit lui-même les variables d'environnement, `PATH`, `ORCADIR` et `XTBEXE` : il n'y a rien à régler à la main, mais les changements ne prennent effet que dans une **fenêtre de commande ouverte après** l'installation. Vérification en une commande, `where.exe orca`. Les modules AUTOCI, s'ils sont nécessaires, se décompressent dans le dossier d'installation et doivent correspondre à la variante installée, série ou parallèle.

Trois pièges spécifiques à Windows. Les calculs se lancent dans un dossier de travail séparé, jamais dans le dossier d'installation. En parallèle, ORCA s'appelle avec son **chemin absolu**, et jamais par `mpirun` : le nombre de cœurs se déclare dans le fichier d'entrée (`!PAL4`, ou un bloc `%pal NProcs ... end` au-delà de huit), et ORCA lance lui-même ses processus. Enfin, un éditeur qui ajoute une extension cachée produit un `input.inp.txt` que rien ne lit.

**Le vrai point de décision sous Windows est le parallélisme** — et il ne se pose que si ORCA est conservé (voir 15.6). ORCA 6 s'appuie sur OpenMPI, et la version attendue est celle qu'indique la page de téléchargement du paquet récupéré (4.1.6 pour la série 6). La documentation est explicite sur Linux et beaucoup plus pauvre sur Windows, où elle finit par suggérer un environnement de type Cygwin. Pour le volume de calculs du projet — des milliers de points de balayage (sections 15.5 et 16) —, faire tourner ORCA sous **WSL2** avec la version Linux est plus simple à maintenir qu'une installation parallèle native Windows, et c'est aussi là que Psi4, xtb et `torsiondrive` s'installent sans friction. Le reste de la chaîne, NAMD et `psfgen`, fonctionne nativement sous Windows : seule la partie quantique gagne à passer par WSL2.

**2. Déclarer le niveau de théorie.** `QuantumSpec` (dans `pepforge.tools.quantum`) porte le logiciel, la méthode, la base, le type de balayage et son pas. Il naît **vide** : `missing_decisions()` liste ce qui reste à fixer et `label()` refuse de produire une étiquette tant que la spécification est incomplète. C'est volontaire, puisque deux entrées de la bibliothèque dérivées à des niveaux différents ne sont pas comparables entre elles.

La spécification quantique n'est **pas** dans `PrepConfig` : une préparation de dynamique n'appelle jamais de quantique. Elle appartient au flux de construction de la bibliothèque (étape 11D), et la valeur employée est recopiée dans chaque profil produit et dans l'entrée de résidu correspondante.

**3. Écrire une implémentation par logiciel.** Le protocole `QuantumTool` est déjà en place avec `optimize`, `scan_dihedral`, `charges` et `version`. Une implémentation est un fichier dans `pepforge/tools/`, qui construit le fichier d'entrée, passe par `ExternalTool` pour que la version soit enregistrée et la sortie analysée, et retourne un `TorsionProfile`. Le reste de PepForge ne change pas d'une ligne.

Deux règles pour ces implémentations. D'abord, **un calcul non convergé lève `ExternalToolError`** et n'entre jamais dans une dérivation : le code de retour d'ORCA ou de Psi4 ne suffit pas, il faut analyser la sortie, et les motifs d'échec se construisent en provoquant délibérément l'erreur sur le logiciel installé (section 17). Ensuite, chaque implémentation passe `tests/contracts.py::assert_quantum_tool_contract`, qui vérifie ce qui ne dépend pas du logiciel : version déclarée, profil non vide, spécification recopiée dans le profil, et profil lisible en relatif.

**4. Produire les profils avec leur provenance.** `TorsionProfile` embarque le fragment, les quatre indices d'atomes de la torsion, les points, la spécification, la version du logiciel, la date et le matériel, et `relative_energies_kcal()` donne les énergies relatives au minimum, seule grandeur comparable au champ de force. Ces profils vivent dans `tests/reference/` et sont ce qui sera montré en soutenance pour justifier les paramètres (section 18.8).

### 15.6 Travailler avec Psi4 seul : oui, et c'est défendable

Question de Firas le 27 novembre 2026 : peut-on se passer d'ORCA, dont l'installation est plus lourde ? **Oui.** C'est un choix légitime, et il a même un avantage méthodologique que l'option à deux logiciels n'a pas.

**Ce que le projet demande réellement au quantique.** Des optimisations de fragments capés de trente à cinquante atomes, des balayages dièdres relaxés unidimensionnels au niveau fixé, éventuellement des charges, et au plus une grille bidimensionnelle sur **un** résidu modèle (section 16). À cette taille et à ce niveau, Psi4 est parfaitement dimensionné : MP2 avec ajustement de densité y est efficace, et c'est exactement le niveau du protocole publié de CGenFF auquel la bibliothèque doit s'aligner.

**L'avantage méthodologique, qui est le vrai argument.** Deux entrées de la bibliothèque dérivées avec deux programmes différents, même au même niveau nominal, ne donnent pas des nombres strictement identiques : les seuils d'intégration, l'ajustement de densité et les critères de convergence diffèrent. **Un seul programme pour toute la bibliothèque, avec sa version enregistrée, est donc plus homogène et plus facile à défendre** qu'une bibliothèque dérivée pour moitié avec chacun. L'instinct de simplifier va ici dans le sens de la rigueur, pas contre elle.

**Ce que l'on perd.** La vitesse brute sur les gros calculs, où ORCA est devant, et l'accès commode à des références de très haute précision du type DLPNO-CCSD(T) si le niveau de théorie devait monter. Aucun des deux n'est nécessaire au périmètre décrit ci-dessus. Si l'un le devenait, ORCA s'ajoute comme **une implémentation de plus** derrière `QuantumTool`, sans toucher au reste : c'est précisément ce que l'interface protège.

**Conséquence pratique immédiate : plus besoin de WSL2.** Psi4 est distribué par conda-forge pour **win-64** en plus de Linux et macOS, et l'installation native Windows est prise en charge depuis sa version 1.4. La recommandation de passer la partie quantique sous WSL2 (section 15.4) portait sur le parallélisme d'ORCA : en abandonnant ORCA, le problème disparaît. NAMD et `psfgen` étant déjà natifs Windows, toute la chaîne tourne alors sur un seul environnement.

```
conda install conda-forge::psi4
```

**Décision retenue : Psi4 seul**, installé par conda-forge, avec sa version enregistrée dans chaque profil produit. ORCA reste documenté en 15.4 et possible plus tard, sans que cela remette en cause une seule dérivation déjà faite, à condition de ne jamais mélanger les deux au sein d'une même bibliothèque.

### 15.5 Automatiser davantage : les briques existantes

Écrire une implémentation par logiciel (section 15.4) n'est pas la seule voie, et ce n'est pas la plus automatique. Quatre briques existent, toutes libres et installables par `pip` ou `conda`, et trois méritent d'être reprises plutôt que réécrites.

**1. Une seule couche de pilotage au lieu d'une implémentation par logiciel : QCEngine.** QCEngine (écosystème MolSSI / QCArchive) expose **une même API Python** pour ORCA, Psi4, PySCF, xtb, Gaussian et d'autres : on décrit la molécule et le calcul voulu, la couche construit le fichier d'entrée, lance le programme et retourne un résultat structuré, provenance comprise. Conséquence pour PepForge : `QuantumTool` reçoit **un seul** adaptateur, celui de QCEngine, au lieu d'un par logiciel, et l'analyse de sortie — la partie fragile et spécifique à chaque version — n'est plus à notre charge. C'est l'option la plus rentable des quatre.

**2. Les balayages de torsion : torsiondrive, et c'est un vrai problème technique, pas du confort.** Un balayage relaxé naïf avance point par point en partant de la géométrie précédente. Il reste alors piégé dans une branche conformationnelle et retourne une courbe qui dépend du sens de parcours : c'est l'**hystérésis**, et elle donne une barrière fausse. `torsiondrive` résout cela par **propagation de front d'onde** : chaque point est relancé depuis toutes les géométries voisines déjà optimisées, et seule la plus basse est conservée. Pour la torsion `N–N`, dont tout le projet dépend, c'est exactement le piège à éviter. La méthode est publiée, et c'est aujourd'hui la façon correcte de produire un profil de torsion de référence.

**3. Du quantique dans le même processus, sans exécutable à trouver : PySCF.** PySCF s'installe par `pip`, s'importe comme un module et n'a ni fichier d'entrée ni fichier de sortie à analyser ; les optimisations de géométrie passent par `geomeTRIC`, qui est aussi l'optimiseur sous contrainte utilisé par `torsiondrive`. Deux avantages pour ce projet : les tests de PepForge peuvent exécuter un vrai calcul quantique minuscule en CI, ce qu'aucun exécutable sous licence ne permet ; et une variante accélérée sur GPU existe, ce qui fait servir deux fois le GPU qui sera de toute façon acheté pour la dynamique (section 16). L'interface Python de Psi4 offre la même propriété.

**4. Une pré-passe semi-empirique : xtb / GFN2.** Les méthodes de liaison fortes semi-empiriques coûtent quelques secondes là où MP2 coûte des minutes. Elles ne remplacent **jamais** les valeurs de référence, mais elles pré-optimisent les géométries, explorent les conformères et éliminent les points inutiles avant la série au niveau fixé. Pratiquement, cela réduit le coût de la série finale sans rien changer à son niveau de théorie, à condition d'écrire noir sur blanc que les valeurs de référence sont toutes au niveau fixé, et jamais semi-empiriques.

**Ce qu'il ne faut pas attendre de ces briques.** Deux limites honnêtes. D'abord, un gestionnaire de file de calculs avec base de données (`QCFractal`) n'a pas sa place ici tant que les séries tiennent sur un poste de bureau : les fichiers de référence de la section 18.8 suffisent, et une base à administrer serait un coût net. Ensuite, les outils d'**ajustement** automatique de paramètres de cet écosystème (`ForceBalance`, chaîne Open Force Field) sont conçus autour de formats de type AMBER/SMIRNOFF, pas de CHARMM : côté quantique ils sont indifférents au champ de force et donc utilisables, mais l'ajustement des paramètres CHARMM36m reste à faire dans PepForge, où la comparaison de profils vit déjà. Enfin, les jeux de données publics de QCArchive contiennent des milliers de balayages de torsion déjà calculés, mais il est très improbable qu'un squelette aza-β³ y figure : à vérifier en quelques minutes, sans compter dessus.

**Recommandation.** Garder ORCA en exécutable externe, parce que sa licence interdit de le distribuer et qu'aucune couche n'y change rien, mais **piloter par QCEngine** plutôt que d'écrire un adaptateur par logiciel, **produire les balayages par torsiondrive**, et **ajouter PySCF** comme implémentation installable par `pip` qui rend les tests exécutables en CI. La pré-passe `xtb` vient ensuite, si le coût de la série le justifie. L'interface `QuantumTool` de la section 15.4 ne change pas : ces briques sont ce qui se trouve derrière elle, et c'est précisément ce que l'interface protège.

**Point de prudence sur le coût.** Des balayages bidimensionnels pour une trentaine d'entrées représenteraient des centaines de milliers de points de calcul. Le périmètre réaliste est : balayages **unidimensionnels** sur les torsions critiques pour chaque entrée, et traitement approfondi sur **un ou deux résidus de référence**. Ce périmètre est à valider avec Firas avant de lancer la moindre série.

---

## 16. Faisabilité des calculs sur un poste de bureau

Firas a posé deux questions le 27 novembre 2026 : combien coûte l'option C de la section 13.3, et l'outil doit-il tourner sur un PC de bureau plutôt que sur un calculateur ? Réponse : oui, à condition de choisir le bon périmètre. Voici les chiffres et la contrainte qui décide.

### 16.1 L'obstacle n'est pas le temps de calcul, c'est le format

Un acide aminé α a **deux** torsions de squelette, `φ` et `ψ` : c'est pourquoi CMAP est une grille bidimensionnelle. Le squelette aza-β³ `-NH-N(R)-CH2-CO-` en a **trois** : autour de `N–N`, de `N–CH2` et de `CH2–C(=O)`. L'équivalent exact de CMAP pour ce squelette serait donc une grille **tridimensionnelle**.

Or le terme CMAP de CHARMM et de NAMD est, sauf erreur à confirmer dans l'étape bibliographique du fil aza, une grille **2D** déclarée sur huit atomes couvrant deux dièdres consécutifs. **Une correction 3D n'est pas exprimable dans ce format.** Deux vérifications à faire avant toute dépense de calcul, et elles sont gratuites :

1. Le motif de huit atomes exigé par CMAP peut-il être déclaré sur la suite d'atomes du squelette aza-β³ ?
2. Si oui, lesquelles des trois torsions le couple ainsi formé décrit-il réellement ?

Si la réponse à la première est non, l'option C tombe, et il reste l'option B. Si elle est oui, l'option C n'est possible qu'en 2D : on retient les deux torsions dominantes dans la grille, et la troisième est traitée par des termes dièdres ordinaires validés par balayage 1D.

### 16.2 Ordres de grandeur

Hypothèses, à calibrer (section 16.4) : résidu aza coiffé à petite chaîne latérale, 25 à 45 atomes ; poste de bureau de 8 à 16 cœurs et 32 à 64 Go ; grille à 15°, soit 24 valeurs par torsion ; « niveau modeste » = DFT en petite base, 2 à 5 min par point optimisé sous contrainte ; « niveau de référence » = DFT dispersive en triple zêta ou RI-MP2, 20 à 60 min par point.

| Travail | Points | Niveau modeste | Niveau de référence |
|---|---|---|---|
| Balayages 1D, 3 torsions, **un** résidu | 72 | 3 à 6 h | 1 à 3 jours |
| Balayages 1D, **toute la bibliothèque** (~30 entrées) | 2 160 | 5 à 8 jours | 1 à 3 mois |
| Grille 2D 24×24, **un** résidu | 576 | 1 à 2 jours | 8 à 24 jours |
| Grille 3D grossière 12×12×12, **un** résidu | 1 728 | 4 à 6 jours | 3 semaines à 2 mois |
| Grille 3D 24×24×24, **un** résidu | 13 824 | 1 à 1,5 mois | 6 mois à 1,5 an |
| Grille 2D × 30 entrées | 17 280 | 1 à 2 mois | 8 mois à 2 ans |

### 16.3 Verdict, et la porte de sortie

| Scénario | Sur un poste de bureau |
|---|---|
| Options A puis B, balayages 1D sur toute la bibliothèque | **faisable** : quelques jours au niveau modeste, et le travail se découpe en lots indépendants qui tournent la nuit |
| Option C en 2D sur **un** résidu de référence | **faisable** : 1 à 2 jours au niveau modeste, jusqu'à trois semaines au niveau de référence |
| Option C en 3D, même grossière | **à écarter** : des semaines pour un seul résidu, et non exprimable dans NAMD |
| Option C sur toute la bibliothèque | **non**, ni sur un poste de bureau ni sur un calculateur modeste |

**La porte de sortie est méthodologiquement propre** : CHARMM n'a pas une grille par résidu, il utilise pour l'essentiel **une seule correction pour tous les acides aminés standards**, avec des cas particuliers pour la glycine et la proline. La même logique s'applique ici : dériver **une** correction de squelette sur un résidu aza modèle à petite chaîne latérale, la transférer aux autres entrées, et valider chaque entrée par des balayages 1D. Le coût total tombe à une grille 2D plus 72 points par résidu, ce qui tient sur un poste de bureau. Cas particuliers à prévoir : l'aza-glycine et tout résidu cyclique.

### 16.4 Calibrer avant d'engager

Les chiffres de 16.2 sont des ordres de grandeur, pas des mesures. Première tâche de l'étape 11D, avant toute série : lancer **une seule** optimisation sous contrainte sur le poste réel, au niveau de théorie envisagé, mesurer le temps de paroi, et multiplier par les nombres de points du tableau. Une heure de travail qui remplace toutes les estimations par un chiffre. Le résultat est consigné dans `tests/reference/` avec la description du matériel.

### 16.5 Le vrai poste de calcul, c'est la dynamique

À ne pas perdre de vue : le quantique n'est pas le gros consommateur du projet. Un peptide de 100 résidus solvaté représente 30 000 à 50 000 atomes. Sur un poste sans accélérateur graphique, la production se compte en quelques nanosecondes par jour ; avec un GPU correct, plutôt en dizaines de nanosecondes par jour. **Pour NAMD, le GPU pèse bien plus lourd que le nombre de cœurs.** Conséquences pratiques :

- l'investissement matériel utile est un GPU, pas un processeur à grand nombre de cœurs ;
- les durées de production et le nombre de répliques (`ProtocolSpec`, étape 9 du plan) doivent être choisis en connaissance de ce débit, et non copiés d'articles utilisant des calculateurs ;
- le rapport consigne le matériel utilisé, pour que les durées de simulation soient interprétables.

### 16.6 Produire les fichiers ne coûte rien ; les défendre coûte

Point à ne pas confondre, parce que la confusion est facile : **les options A, B et C de la section 13.3 ne concernent que la correction CMAP, pas la production des fichiers.** Aucune n'est un prérequis pour obtenir un `.str`, un `.prm`, un `.psf` et un `.pdb` de pseudopeptide.

| Travail | Ce qu'il demande | Coût sur un poste de bureau |
|---|---|---|
| `.str` du fragment aza par CGenFF | CGenFF attribue les types et les charges **par analogie**, avec des pénalités. Ce n'est pas un calcul quantique | secondes à une minute |
| Transformation en résidu de chaîne : `RESI`, patches `PRES`, charges de jonction entières | travail de **développement** (étape 11C), pas de calcul | aucun calcul |
| `.prm` | CHARMM36m fourni, plus les ajouts dérivés | aucun calcul |
| `.psf` et `.pdb` par `psfgen` | assemblage du système | secondes |
| **Valider** que ces fichiers donnent la bonne énergétique de squelette | balayages 1D, options A puis B | heures à jours par résidu |
| Option C, grille CMAP | voir 16.2 et 16.3 | 1 à 2 jours pour un résidu modèle, au niveau modeste |

Autrement dit : produire les fichiers est immédiat, et c'est précisément le danger. **Un `.psf` incorrect s'assemble et se simule sans rien signaler** ; la trajectoire aura l'air normale alors que les préférences conformationnelles du squelette modifié seront fausses. C'est pour cela que la validation existe, et que les contrôles de la section 13 énumèrent les termes depuis le `.psf` réellement produit : un paramètre manquant doit lever une erreur, jamais être omis en silence.

Ordre de travail recommandé, qui va du moins coûteux au plus coûteux :

1. produire les fichiers (minutes) ;
2. faire passer les contrôles d'intégration de la section 13.1 (aucun calcul) ;
3. lancer les balayages 1D et mesurer l'écart (heures à jours) ;
4. **au vu des chiffres seulement**, décider si l'option B suffit ou si l'option C se justifie.

---

## 17. Politique d'absence d'erreur silencieuse

Exigence posée par Firas le 27 novembre 2026 : **aucune erreur silencieuse, ni dans le `.psf` ni dans les autres fichiers**. C'est la plus importante du projet, parce qu'une erreur silencieuse ne se voit pas à la lecture d'une trajectoire : elle se voit à la soutenance.

### 17.1 Le principe : ne jamais croire l'outil sur parole

Un code de retour nul ne prouve rien, et un avertissement affiché dans un flot de sortie de plusieurs centaines de lignes n'est pas un signalement. La seule méthode fiable est l'**inventaire attendu contre fichier produit** :

1. Avant d'écrire quoi que ce soit, PepForge construit, à partir du `Peptide` et du `ParameterSet`, un **inventaire attendu** : nombre d'atomes par résidu, résidus dans l'ordre, liaisons du squelette, extrémités, ponts disulfure, charge nette, nombre de termes CMAP attendus.
2. Après écriture, PepForge **relit le fichier réellement produit** et en extrait le même inventaire.
3. Les deux inventaires sont comparés terme à terme. Toute différence lève `QualityCheckError`.

Cette méthode ne dépend ni de la version de `psfgen`, ni de la façon dont il formule ses avertissements. C'est ce qui la rend robuste.

**Ce mécanisme est implémenté** dans [`pepforge.checks.inventory`](../src/pepforge/checks/inventory.py), avec ses tests dans [`tests/test_inventory.py`](../tests/test_inventory.py) : un test par piège du catalogue 17.2. Deux décisions de conception y portent toute la valeur :

1. **Un champ non vérifiable est un échec, jamais un succès.** Si l'une des deux faces ignore une valeur, le résultat est `FAIL` avec la mention « non vérifiable ». « Je n'ai pas pu contrôler » ne doit jamais se lire comme « tout va bien » : c'est exactement par là qu'une erreur silencieuse s'installe.
2. **La comparaison ne s'arrête pas au premier écart** : un résultat par champ, pour qu'un rapport montre tous les problèmes d'un coup.

L'extraction de l'inventaire d'un fichier produit est propre à chaque moteur et vit dans `engines` ; la comparaison, elle, est pure, et se teste donc sans aucun outil externe installé.

**Depuis l'étape 1**, les ponts disulfure de l'inventaire attendu viennent des liaisons déclarées sur le `Peptide`. Une déclaration qui n'a pas de champ dans l'inventaire (cyclisation tête-queue, modification de chaîne latérale) fait **refuser** la construction de l'inventaire attendu (`UnsupportedFeatureError`) : un inventaire qui l'omettrait laisserait écrire un fichier sans elle, et la comparaison ne pourrait pas s'en apercevoir. Le champ correspondant s'ajoutera avec la prise en charge de la déclaration, pas avant.

### 17.2 Catalogue des défaillances silencieuses

Ce catalogue est à compléter au fil des outils réellement utilisés. Chaque ligne devient un contrôle nommé dans `checks`, avec son test `PASS` et son test `FAIL`.

| Étape | Ce qui se passe sans contrôle | Contrôle qui l'attrape |
|---|---|---|
| Lecture PDB | conformations alternatives non résolues : un atome sur deux est retenu arbitrairement | inventaire : aucun `alt_loc` restant |
| Lecture PDB | résidus manquants dans la structure, donc **rupture de chaîne** au milieu du peptide | continuité de la numérotation et des liaisons du squelette |
| `psfgen` | un nom de résidu absent de la topologie : le résidu peut être **ignoré**, et le `.psf` sort amputé | comparaison du nombre de résidus et d'atomes attendus contre produits |
| `psfgen` | un atome du PDB absent de la topologie, ou l'inverse : écarté ou régénéré | comparaison atome par atome, par nom et par résidu |
| `psfgen` | coordonnées manquantes devinées automatiquement, parfois mal placées | aucun atome à une position dupliquée ou aberrante ; distances de liaison dans une plage sourcée |
| `psfgen` | patch non appliqué ou appliqué au mauvais résidu : extrémités fausses, jonction absente | les extrémités produites correspondent à celles du `Peptide` ; chaque jonction attendue porte sa liaison |
| `psfgen` | pont disulfure non déclaré : deux cystéines libres au lieu d'un pont | inventaire des ponts attendus contre liaisons `S–S` présentes |
| `psfgen` | **terme CMAP absent** pour un résidu standard, parce qu'une référence traversant la jonction ne s'est pas résolue | comptage des termes CMAP du `.psf` contre le nombre de résidus standards attendus (section 13.3) |
| Paramètres | un dièdre sans paramètre propre est capté par une **entrée générique à jokers** : un paramètre s'applique, mais ce n'est pas le bon | détection et comptage des termes résolus par joker, listés dans le rapport |
| Paramètres | paramètres 1-4 explicites non déclarés : les valeurs non liées ordinaires s'appliquent à la place | vérification que chaque type concerné déclare `eps14` et `rmin14` là où CHARMM36m le fait |
| Paramètres | termes `NBFIX` non appliqués parce que les types diffèrent | liste des `NBFIX` attendus contre ceux effectivement actifs |
| Charges | redistribution de jonction approximative : charge nette non entière de quelques centièmes | charge nette relue dans le `.psf`, comparée à la charge théorique, tolérance stricte |
| Solvatation | chevauchement atomique, ou marge de boîte insuffisante | distance minimale entre atomes non liés ; marge mesurée sur le système produit |
| Protonation | état choisi par un outil tiers sans que ce soit demandé | chaque état de protonation du système produit est comparé à celui décidé à l'étape 2 |
| Conversion entre moteurs | perte des termes CMAP, des paramètres 1-4 ou altération des charges | nombre d'atomes, charge nette **et énergie en point unique** comparés avant et après |
| Métal | centre métallique dissocié du chélateur en cours de simulation sans que rien ne le signale | géométrie de coordination vérifiée après construction, puis contrôlée en fin d'équilibration |

### 17.3 Les avertissements sont des erreurs

Règle par défaut : **toute ligne d'avertissement d'un outil externe fait échouer l'étape.** Une exception ne s'obtient qu'à trois conditions, cumulatives :

1. le motif exact est inscrit dans une liste d'exceptions versionnée, dans les données du projet ;
2. la raison est écrite, avec la version de l'outil concernée ;
3. l'exception est reportée dans le rapport de chaque préparation où elle s'applique, et non seulement dans le code.

La liste des motifs d'échec de `ExternalTool` ne peut pas être devinée : elle doit être **construite et vérifiée contre les outils réellement installés**, en provoquant délibérément chaque erreur pour observer ce que l'outil affiche. C'est un travail d'étape, pas une liste écrite d'avance.

### 17.4 Contrôles de dernier recours

Même avec tout ce qui précède, trois contrôles numériques attrapent ce que l'inventaire ne voit pas :

- **Charge nette relue** dans les fichiers produits, jamais recopiée de l'entrée.
- **Énergie en point unique** calculée sur le système produit : une valeur aberrante révèle un paramètre absurde ou un chevauchement, là où les comptages sont tous corrects.
- **Minimisation courte qui doit converger** : un système mal paramétré diverge ou stagne.

Ces trois valeurs sont consignées dans le rapport, avec le matériel et les versions d'outils, pour qu'une préparation reste comparable à une autre.

### 17.5 Ce que cela impose au code

- `EngineFiles` porte `net_charge` et `atom_count` **relus dans les fichiers produits** (section 4) : c'est déjà le contrat, et cette section explique pourquoi.
- Chaque contrôle du catalogue 17.2 est une fonction pure de `checks`, avec son test `PASS` et son test `FAIL` (règle 7 du projet).
- `ExternalTool` analyse la sortie et non le seul code de retour : c'est déjà son contrat, la liste de motifs reste à construire.
- Un échec arrête le pipeline par `QualityCheckError` et produit un **rapport partiel lisible** : savoir où ça s'est arrêté fait partie de l'absence d'erreur silencieuse.

---

## 18. Comment vérifier et valider une sortie

La section 7 décrit comment on teste le **code**. Celle-ci décrit comment on contrôle une **sortie scientifique** : les fichiers d'un pseudopeptide donné. Question posée par Firas le 27 novembre 2026.

D'abord, deux mots qu'il ne faut pas confondre, parce qu'ils coûtent des choses très différentes :

- **Vérifier** : le fichier dit-il ce que je voulais dire ? C'est une question de conformité. Gratuit, automatisable, et c'est ce qui attrape les erreurs silencieuses.
- **Valider** : le modèle décrit-il correctement la molécule ? C'est une question de pertinence physique. Coûteux, et cela exige une référence extérieure, quantique ou expérimentale.

Un fichier peut être parfaitement vérifié et complètement invalide. C'est le cas le plus dangereux, et c'est précisément celui du résidu aza-β³ paramétré par analogie.

### 18.1 L'échelle, du moins coûteux au plus coûteux

| Niveau | Ce qu'on contrôle | Coût | Ce qu'il attrape |
|---|---|---|---|
| 0 | conformité des fichiers | minutes | erreurs silencieuses de construction |
| 1 | cohérence numérique | minutes | paramètres absurdes, chevauchements |
| 2 | cohérence structurale | heures | géométrie fausse du résidu modifié |
| 3 | cohérence énergétique | jours | énergétique de squelette fausse |
| 4 | accord expérimental | la thèse | le modèle décrit-il le vrai peptide |

Un niveau ne se saute pas : chacun attrape des défauts que le suivant masquerait.

### 18.2 Niveau 0 — conformité des fichiers

À faire pour **chaque** sortie, sans exception. C'est ce que `pepforge.checks.inventory` automatise, et que l'on peut aussi contrôler à la main.

1. **Charge nette** : sommer les charges partielles du `.psf` et comparer à la charge théorique attendue au pH retenu. Dans VMD, la somme pondérée par la charge d'une sélection donne cette valeur directement. Tolérance stricte : un écart de quelques centièmes trahit une redistribution de jonction approximative.
2. **Nombre d'atomes et de résidus**, comparés à ce qui était voulu, résidu par résidu. Un résidu ignoré par l'outil se voit ici, et nulle part ailleurs.
3. **Le résidu modifié existe** sous le nom attendu, avec son nombre d'atomes attendu.
4. **Les liaisons de jonction existent.** Visuellement, dans VMD, une liaison de squelette absente apparaît comme une rupture dans la chaîne. Par script, on vérifie que les deux atomes de jonction figurent bien dans la liste des liaisons l'un de l'autre.
5. **Les extrémités** produites correspondent à celles demandées.
6. **Les termes croisés** : la section du `.psf` qui porte les termes CMAP (`NCRTERM`) doit en compter autant qu'il y a de résidus standards attendus. Un terme manquant ne provoque aucune erreur à l'exécution (section 13.3).
7. **Aucune coordonnée aberrante** : pas d'atome à l'origine, pas deux atomes au même endroit. C'est la trace d'une coordonnée devinée par l'outil.

### 18.3 Niveau 1 — cohérence numérique

8. **Charger le système dans NAMD pour zéro pas** et lire les énergies initiales. C'est le test le plus rentable du lot : un paramètre lié absent fait échouer NAMD bruyamment, et une énergie absurde ou non finie saute aux yeux. Quelques secondes.
9. **Minimisation courte** : l'énergie doit décroître régulièrement puis se stabiliser. Une divergence, une stagnation immédiate ou un atome qui part à l'infini signalent un paramètre faux ou un chevauchement.
10. **Comparer la géométrie avant et après minimisation.** Si le déplacement est concentré sur le résidu modifié, ses paramètres sont en cause ; si tout bouge, c'est la structure de départ.

### 18.4 Niveau 2 — cohérence structurale

11. **Dynamique courte**, puis contrôle que la chimie du résidu modifié tient : longueurs de liaison stables, liaison `N–N` qui ne s'allonge pas anormalement, amide `ω` restant proche de la planéité, chiralité des autres résidus préservée.
12. **Inspection visuelle** du voisinage modifié dans VMD ou PyMOL. À ne pas négliger : l'œil repère des choses qu'aucun contrôle ne cherche.
13. **Comparer la géométrie du résidu aza** obtenue en champ de force à la géométrie **optimisée en quantique** du fragment coiffé. C'est la première confrontation à une référence extérieure.

### 18.5 Niveau 3 — cohérence énergétique

14. **Profils de torsion** champ de force contre quantique sur les torsions critiques, dont `N–N` (étapes 11D et 13.4).
15. **Répliques indépendantes** : deux simulations aux vitesses initiales différentes doivent converger vers les mêmes populations conformationnelles. Si elles divergent, l'échantillonnage est insuffisant, et aucune conclusion n'est permise — ce n'est pas nécessairement une erreur de paramètres, mais cela interdit de conclure.
16. **Peptide témoin, et c'est le contrôle le plus rentable de toute la liste.** Simuler le peptide **non modifié**, en CHARMM36m pur, dans exactement les mêmes conditions. S'il se comporte conformément à ce qui est publié pour lui et que la version aza part dans tous les sens, le problème est dans les paramètres dérivés, pas dans le protocole. Et inversement : si le témoin se comporte déjà mal, inutile d'accuser le résidu modifié. Ce témoin coûte une simulation et sépare proprement les deux sources d'erreur.

### 18.6 Niveau 4 — accord expérimental

17. **Confrontation aux données RMN** du peptide, que le projet de thèse prévoit de produire (CARA, ARIA) : distances issues des NOE, couplages, déplacements chimiques secondaires. Pour un peptide, c'est la validation décisive, et elle sera disponible en interne.
18. Pour le complexe, confrontation du comportement simulé aux mesures d'affinité et de spécificité.

### 18.7 Deux règles de méthode

- **Les seuils se fixent avant de regarder les résultats**, et s'écrivent. Un seuil choisi après coup transforme la validation en justification.
- **Les valeurs de référence vivent dans `tests/reference/`**, avec leur provenance et le matériel utilisé, afin qu'une sortie reste comparable à une autre et qu'un contrôle puisse être rejoué des mois plus tard.

### 18.8 Exemple complet de référence quantique

Une « référence quantique » n'est pas un nombre isolé : c'est une **courbe d'énergie** calculée sur un fragment défini, à un niveau de théorie fixé, que le champ de force doit reproduire. Voici l'exemple de bout en bout pour la torsion la plus critique du projet, la liaison `N–N` du squelette aza-β³.

**1. Le fragment, et non le peptide.** On ne fait pas de quantique sur un peptide de 20 résidus. On construit le plus petit fragment qui porte la torsion avec son environnement électronique réel : le résidu aza-β³ **capé** des deux côtés par des groupes acétyle et *N*-méthylamide, soit `CH₃-CO-NH-N(R)-CH₂-CO-NH-CH₃`. Les capes imitent les liaisons peptidiques voisines ; sans elles, la torsion mesurée serait celle d'une amine terminale, pas celle d'un squelette.

**2. Le balayage.** Balayage **relaxé** de la torsion `C(=O)–N–N–CH₂` : on contraint l'angle à une valeur, on optimise tout le reste, on relève l'énergie, et on recommence. 36 points tous les 10°, soit un tour complet. Relaxé et non rigide : un balayage rigide attribue à la torsion le coût de recouvrements stériques qu'une molécule réelle évite en se réarrangeant.

**3. Le niveau de théorie, fixé une fois pour toute la bibliothèque.** Deux entrées dérivées à des niveaux différents ne sont pas comparables (section 15). Le protocole publié de CGenFF utilise des balayages au niveau MP2 avec une base de type 6-31G(d) : s'aligner sur lui est l'option la plus défendable, puisque les paramètres dérivés doivent cohabiter avec les siens. **À confirmer sur la publication de référence avant de lancer la première série**, et à inscrire dans la configuration. Une alternative moderne, moins coûteuse à précision comparable sur des torsions, est une fonctionnelle à dispersion corrigée avec une base triple-zêta ; le choix est à arrêter avec Firas, pas par défaut.

**4. L'entrée, concrètement.** Avec ORCA, le balayage ci-dessus s'écrit ainsi (les indices d'atomes sont ceux du fichier de géométrie, et sont écrits par PepForge, jamais à la main) :

```
! MP2 6-31G(d) Opt TightSCF
%geom Scan
  D 1 2 3 4 = 0.0, 350.0, 36
  end
end
* xyzfile 0 1 fragment_aza.xyz
```

Avec Psi4, le même balayage se pilote depuis Python, ce qui est la raison de le garder en second outil derrière la même interface `QuantumTool` : la boucle sur les 36 points, l'extraction des énergies et l'écriture du fichier de référence tiennent dans le même script que le reste de PepForge.

**5. Ce à quoi on compare.** Le même balayage est refait **avec le champ de force**, sur le même fragment, dans NAMD : même torsion contrainte, mêmes 36 valeurs, minimisation à chaque point. On obtient deux courbes. Pour l'**ajustement** du paramètre dièdre, ce n'est pas la courbe quantique qu'on reproduit directement mais la **différence** entre la courbe quantique et la courbe du champ de force calculée *sans* le terme dièdre en cours d'ajustement : le reste des termes (angles, non liées, 1-4) contribue déjà au profil, et l'ignorer reviendrait à compter deux fois la même énergie.

**6. Le critère de réussite, écrit avant de regarder.** Trois conditions, à fixer avec Firas puis à inscrire dans `tests/reference/` :

- les **minima** des deux courbes au même endroit, à une tolérance angulaire près ;
- l'**ordre** des minima préservé, c'est-à-dire que le conformère le plus stable en quantique reste le plus stable avec le champ de force ;
- l'**écart quadratique moyen** sur la fenêtre énergétique réellement peuplée à température ambiante sous un seuil donné, exprimé en kcal/mol.

Le troisième critère ne se juge que sur la fenêtre peuplée : reproduire le sommet d'une barrière à 15 kcal/mol n'a aucune conséquence sur une trajectoire, alors qu'une erreur de 1 kcal/mol entre deux minima proches change les populations.

**7. Ce que le fichier de référence contient.** Pour chaque entrée de la bibliothèque, un fichier versionné portant : la géométrie du fragment, la torsion balayée par ses indices d'atomes, le niveau de théorie, le logiciel et sa version, les 36 énergies, la date, et le matériel. C'est ce qui rend le contrôle rejouable des mois plus tard, et c'est ce qui sera montré en soutenance pour justifier les paramètres.

### 18.9 Choisir le niveau de validation dans l'outil

Le niveau est un **choix de l'utilisateur**, porté par `PrepConfig.validation_level` (`pepforge.config.ValidationLevel`), parce que le coût varie d'un facteur énorme d'un palier à l'autre. Trois propriétés sont garanties par le code (`pepforge.checks.levels`, testé par `tests/test_levels.py`) :

1. **Les niveaux sont cumulatifs.** Demander le niveau 3 exige aussi les niveaux 0, 1 et 2. Un contrôle énergétique sur un fichier dont la conformité n'a pas été vérifiée ne veut rien dire.
2. **Un contrôle demandé et non exécuté est une erreur.** `audit()` confronte le catalogue des contrôles du niveau demandé aux résultats réellement produits, et rapporte un `FAIL` pour chaque contrôle manquant, en distinguant deux cas : « pas encore implémenté dans PepForge », avec l'étape du plan qui le fournira, et « demandé mais non exécuté par cette préparation ». Sans ce mécanisme, un niveau de validation serait une promesse : l'utilisateur demanderait le niveau 3, l'outil ferait ce qu'il sait faire, et le rapport ne dirait pas que les profils de torsion n'ont jamais été calculés.
3. **Un niveau supérieur ne masque jamais un échec inférieur.** Les résultats des paliers bas sont repris tels quels dans l'audit.

`unavailable(level)` permet d'avertir **avant** de lancer une préparation longue, et non de le découvrir dans le rapport final. État actuel du catalogue : le niveau 0 est entièrement outillé sauf le contrôle de géométrie (étape 4) ; les niveaux 1 à 3 sont déclarés et refusés tant que les étapes 7, 9, 11D et 14 ne sont pas écrites ; le niveau 4 est hors de l'outil par nature, la confrontation aux données RMN restant manuelle.

Jusqu'où est-ce faisable sur un poste de bureau ? **Le niveau 3 est atteignable**, à deux conditions : le palier 3 n'est exécuté que sur **un ou deux peptides de référence**, pas sur chaque préparation, et le peptide témoin de la section 18.5 est simulé une seule fois par protocole, pas une fois par peptide. Le niveau 0 tourne en quelques secondes et peut rester actif en permanence ; les niveaux 1 et 2 coûtent de quelques minutes à quelques heures et se justifient à chaque nouvelle entrée de la bibliothèque ; le niveau 3 coûte des jours, et c'est le palier qui rend les paramètres défendables en soutenance.

---

## 19. Distribution : un seul code, deux exécutables

Question de Firas le 27 novembre 2026 : faut-il une version exécutable qui tourne sous Linux **et** Windows, ou deux versions, une par système ?

### 19.1 La réponse en une ligne

**Un seul code source, jamais deux.** PepForge est du Python pur : le même dépôt tourne sous Linux, Windows et macOS sans la moindre variante de code. En revanche, **un exécutable autonome est propre à un système par nature** : aucun binaire unique ne s'exécute sur les deux. Donc, si exécutable autonome il y a, ce sont deux **constructions** du même code, produites automatiquement par la CI sur chaque système, et non deux versions à maintenir.

Maintenir deux versions du code serait la pire option : deux fois la correction de chaque bogue, et la garantie qu'elles divergent.

### 19.2 Le vrai obstacle n'est pas PepForge, ce sont les outils externes

C'est le point qu'il faut voir avant de choisir un format de distribution. PepForge ne fait pas le calcul : il pilote `psfgen`, NAMD, Psi4 et, plus tard, `gmx` ou les outils AMBER. Ces logiciels ont leurs propres installateurs, leurs propres licences, et une disponibilité qui varie selon le système. PepForge **ne peut pas les embarquer** : les licences l'interdisent pour certains, et le volume l'interdirait pour les autres.

Conséquence directe : **un `.exe` ne dispense pas d'installer les outils externes.** Un exécutable autonome rendrait seulement Python invisible, ce qui est le problème le plus petit de la chaîne. C'est pourquoi la distribution ne commence pas par l'empaquetage mais par la **détection** : `ExternalTool` cherche chaque exécutable sur le `PATH`, enregistre sa version et échoue clairement s'il manque, au lieu de se taire.

### 19.3 Les quatre voies, et laquelle choisir

| Voie | Ce que l'utilisateur fait | Pour | Contre |
|---|---|---|---|
| **`pipx install pepforge`** | une commande, puis `pepforge ...` | la plus simple à publier et à mettre à jour, environnement isolé, fonctionne identiquement sur les deux systèmes | demande Python installé |
| **Paquet conda** | `conda install pepforge` dans un environnement | peut **déclarer Psi4 et xtb comme dépendances** et les installer avec, ce qu'aucune autre voie ne sait faire | demande conda ; ne couvre ni NAMD ni VMD, absents des canaux usuels |
| **Exécutable autonome** (PyInstaller ou Nuitka) | télécharger un fichier et le lancer | ne demande pas Python ; utile pour confier l'outil à un collègue | un binaire **par système**, construit sur ce système ; démarrage plus lent ; masque l'environnement Python, donc diagnostic plus difficile ; ne dispense pas des outils externes |
| **Image conteneur** | `docker run ...` | reproductibilité maximale, utile pour qu'un rapporteur rejoue une préparation des années plus tard | sous Windows, s'appuie sur WSL2 ; les licences de NAMD et VMD interdisent de redistribuer une image qui les contient, donc l'image se construit chez l'utilisateur avec ses propres archives |

**Recommandation pour ce projet** : `pipx` comme voie normale, **paquet conda** ensuite parce que c'est la seule qui installe Psi4 en même temps, et exécutable autonome **seulement si** le besoin apparaît de confier l'outil à quelqu'un qui n'installera pas Python. Le conteneur vient en dernier, et sous forme de recette plutôt que d'image publiée, pour la reproductibilité de la thèse.

### 19.4 Ce que la portabilité exige dans le code

Elle s'obtient par des règles simples, appliquées dès maintenant plutôt que réparées à la fin :

- **Chemins** : `pathlib.Path` partout, jamais de concaténation de chaînes ni de séparateur écrit à la main.
- **Appels externes** : `subprocess` avec une liste d'arguments, jamais `shell=True` ; découverte par `shutil.which`, jamais un chemin codé en dur.
- **Encodage explicite** : `encoding="utf-8"` à chaque lecture et écriture, parce que la valeur par défaut de Windows n'est pas UTF-8.
- **Fins de ligne** : les fichiers destinés à un outil externe (`.pdb`, `.inp`, scripts NAMD) s'écrivent en `\n`, donc avec `newline="\n"` explicite. Un fichier en `CRLF` passe souvent, mais pas partout, et l'échec est précisément du genre silencieux que la section 17 interdit.
- **Casse des noms de fichiers** : jamais deux fichiers ne différant que par la casse, le système de fichiers de Windows ne les distinguant pas.
- **Chemins longs** : garder les dossiers de travail courts, la limite historique de Windows frappant les arborescences profondes.

### 19.5 Comment on s'assure que c'est vrai

**La CI tourne sur les deux systèmes** (`.github/workflows/tests.yml`) : Linux sur les trois versions de Python prises en charge, Windows sur la plus ancienne et la plus récente, là où les écarts de portabilité se voient. La portabilité est ainsi **vérifiée à chaque poussée**, et non supposée. C'est la même logique que la matrice de dépendances : une règle qui n'est pas exécutée n'est pas une règle.

L'empaquetage proprement dit est une étape tardive du plan (étape 15), après que la chaîne fonctionne de bout en bout : empaqueter un outil incomplet n'apporte rien, alors que tester sur les deux systèmes dès maintenant évite d'accumuler des dettes de portabilité.

---

## 20. Entrées, traitement par lots et sorties proposées

Questions de Firas le 27 novembre 2026 : quelles entrées accepter, y compris pour traiter un lot de peptides d'un coup, et quelles sorties proposer au-delà de `.prm`, `.str`, `.psf` et `.pdb` ?

### 20.1 Ce que chaque format d'entrée porte réellement

Le point à retenir : ces formats ne sont pas des syntaxes interchangeables, ils portent des **informations différentes**. Choisir le mauvais format fait perdre de l'information en silence.

| Format | Porte | Ne porte pas | Usage juste dans ce projet |
|---|---|---|---|
| **Séquence / notation** (`Ac-RGDS-NH2`) | enchaînement des résidus, extrémités, modifications déclarées | aucune coordonnée | la voie principale de conception, et **la seule qui passe à l'échelle d'un lot** |
| **PDB** | coordonnées, noms de résidus, chaînes | ordres de liaison ; limité à un identifiant de chaîne d'un caractère et à 99 999 atomes | récepteur, poses de docking, modèles RMN |
| **mmCIF** | tout ce que porte le PDB, sans ses limites de taille, plus les métadonnées et les modèles multiples | ordres de liaison | **à préférer au PDB** pour tout ce qui vient de la PDB ou dépasse une chaîne simple |
| **SDF / MOL** | liaisons avec leurs ordres, hydrogènes explicites, plusieurs molécules par fichier | notion de résidu et de chaîne | fragments : chélateur, résidu nouveau, soumission à CGenFF. **Pas** un peptide de 100 résidus |
| **MOL2** | liaisons, **types atomiques et charges partielles** | — | utile, mais ses types atomiques sont ceux du programme qui l'a écrit, **jamais des types CHARMM** : les adopter serait exactement l'erreur silencieuse de la section 17 |
| **SMILES** | connectivité, stéréochimie si elle est écrite | aucune coordonnée, aucune conformation | définir compactement une chaîne latérale ou un chélateur |

Deux pièges à traiter explicitement dans `io`. D'abord, un fichier **multi-modèles** (les structures RMN à venir en produiront) : le lecteur exige que le modèle soit choisi, et ne prend jamais le premier en silence. Ensuite, les **conformations alternatives** : une occupation partielle non résolue est au catalogue des pièges silencieux (section 17.2).

### 20.2 Deux entrées supplémentaires que je propose

**1. Un manifeste, pour le lot.** C'est la vraie réponse au traitement par lots. Un dossier de deux cents fichiers PDB n'est pas la description d'un lot, parce que les substitutions, les extrémités et le pH ne sont pas dans ces fichiers. Un manifeste CSV ou YAML porte une ligne par peptide :

```
id,sequence,notation,substitutions,n_terminus,c_terminus,ph,goal
P001,Ac-RGDS-NH2,one_letter,,acetyl,amide,7.4,free_peptide
P002,Ac-RGDS-NH2,one_letter,3:aza,acetyl,amide,7.4,free_peptide
P003,Ac-Arg-Gly-Asp-Ser-NH2,three_letter,"1:aza,3:aza",acetyl,amide,7.4,free_peptide
```

La colonne `substitutions` est ce qui fait de PepForge autre chose qu'un convertisseur de fichiers : elle décrit le balayage de positions de l'étape 11E. La colonne `notation` est obligatoire : le mode de lecture d'une séquence n'est jamais deviné (section 5.7).

**2. La notation HELM**, en second format textuel. Elle a été conçue pour les biopolymères à **monomères non naturels**, c'est-à-dire exactement le cas de ce projet, et elle fait l'aller-retour sans perte. Elle vaut d'être lue en entrée, même si la notation interne reste celle de la séquence.

Et la configuration YAML de PepForge est elle-même une entrée : le rapport la régénère (règle 8), donc rejouer une préparation à l'identique ne coûte rien.

### 20.3 Comment un lot doit se comporter

```
pepforge prepare --manifest peptides.csv --outdir runs/
```

Quatre propriétés, et la première est la plus importante :

1. **Isolement** : l'échec d'un peptide n'arrête jamais le lot. Sinon un lot de cinquante s'arrête au septième et la nuit de calcul est perdue.
2. **Un dossier par peptide**, nommé par son `id`, contenant l'intégralité de ses sorties et de son rapport.
3. **Un récapitulatif** `summary.csv` et `summary.md` : par peptide, l'état (réussi, échoué, **contrôles non vérifiés**), la charge nette, le nombre d'atomes, le nombre de substitutions, et les identifiants des contrôles en échec. Les trois états sont distincts : « non vérifié » n'est pas « réussi » (section 18.9).
4. **Code de retour non nul si un seul peptide a échoué**, pour que le lot s'enchaîne dans un script sans qu'on ait à lire la sortie.

### 20.4 Les sorties que je propose, au-delà de `.prm`, `.str`, `.psf` et `.pdb`

**Ce qu'il faut pour lancer réellement le calcul.**

1. **Les scripts NAMD, un par phase** : minimisation, chauffage, équilibration NVT puis NPT avec relâchement progressif des contraintes, production. Avec le **fichier de référence des contraintes** qu'ils utilisent, sans lequel ils ne tournent pas.
2. **Le `.crd` CHARMM** à côté du `.pdb` : il porte plus de précision et certains chemins de `psfgen` et de VMD le préfèrent.
3. **Le système sec et le système solvaté comme deux jeux distincts**, pour pouvoir revenir en arrière sans refaire la paramétrisation.

**Ce qu'il faut pour pouvoir y faire confiance — la partie sur laquelle j'insisterais.**

4. **Le rapport de vérification en deux formes** : lisible (Markdown) et exploitable (JSON), une entrée par contrôle avec son niveau, son état, son message et les deux inventaires comparés. Le JSON se **compare d'une exécution à l'autre**, et c'est ainsi qu'on voit une régression apparaître.
5. **Le journal** des modifications : quoi, pourquoi, avec quels paramètres (`Journal`).
6. **La configuration YAML régénérée**, suffisante à elle seule pour rejouer la préparation.
7. **Les empreintes** de chaque fichier d'entrée et de sortie, et la **version de chaque outil externe**. Sans cela, « la même préparation » ne veut rien dire six mois plus tard.
8. **La table de correspondance des numérotations** : numéro interne du résidu contre numéro dans les fichiers produits. `psfgen` et le format PDB renumérotent ; une correspondance fausse ruine **silencieusement** toutes les analyses ultérieures. Cette sortie coûte presque rien et évite une erreur particulièrement vicieuse.
9. **L'inventaire des termes extrait du `.psf` réellement produit** : liaisons, angles, dièdres, impropres et termes CMAP traversant chaque jonction aza (section 13). C'est déjà une exigence de vérification ; ce doit être un fichier, pas seulement un contrôle, parce que c'est la pièce qui prouve l'intégration.

**Ce qui sert directement la thèse.**

10. **Le fichier des sources** (`.bib` ou liste) : champ de force et sa version, modèle d'eau, chaque entrée de paramètre avec sa source et son niveau de confiance, chaque version d'outil, le niveau de théorie quantique. Il alimente le chapitre méthodes sans recopie manuelle.
11. **Un brouillon de paragraphe « méthodes »** en français, engendré depuis la configuration et le journal. À relire, jamais à publier tel quel, mais il évite les oublis.
12. **Un script VMD** qui charge le système et met en évidence les résidus aza et leurs jonctions, pour l'inspection visuelle du niveau 2 (section 18.4).

**Plus tard.** Le lot de conversion vers GROMACS puis AMBER, comme étape explicite et vérifiée (étapes 13A et 13B), jamais comme sortie implicite.

### 20.5 Quand ces sorties arrivent

Les sorties 1 à 3 viennent avec les étapes 7 à 9. La table de correspondance (8) et l'inventaire des termes (9) viennent avec l'étape 7, parce que c'est là que le `.psf` naît. Le rapport, le journal, la configuration et les empreintes (4 à 7) viennent avec l'étape 10, qui est précisément l'étape du rapport reproductible. Le fichier des sources et le brouillon de méthodes (10, 11) viennent avec l'étape 11E, quand la bibliothèque est complète. Le manifeste et le mode lot s'ajoutent à l'étape 11E également, puisque c'est le balayage de positions qui les rend nécessaires.

---

## 21. Interface utilisateur

Exigence de Firas le 27 novembre 2026 : une interface simple, intuitive, utilisable sans compétence en programmation, au design professionnel, moderne, épuré et cohérent, avec une navigation claire, des formulaires simples, des messages d'erreur explicites, un affichage lisible des résultats et une identité visuelle sobre.

### 21.1 Le principe qui conditionne tout le reste

**L'interface est un client du même pipeline que la ligne de commande, jamais un second chemin de code.** Ni `ui` ni `cli` ne contiennent de logique scientifique : les deux lisent des choix, appellent `pipeline`, et affichent ce qu'il retourne. Une interface qui appellerait `psfgen` directement, ou qui court-circuiterait les contrôles de la section 17 pour « aller plus vite », produirait exactement les erreurs silencieuses que tout le reste de cette architecture s'emploie à rendre impossibles.

Conséquence pratique : **toute action dans l'interface se traduit par la même configuration YAML** que la ligne de commande consomme, et l'interface la montre. La reproductibilité (règle 8) survit ainsi au passage à l'interface graphique, au lieu d'être perdue dans des clics non tracés.

**Dépendances** : un module `ui` au sommet de la matrice, autorisé à importer `pipeline`, `config`, `checks`, `report`, `io` et `core`, et rien d'autre. Sa ligne est ajoutée à la matrice de la section 3 et à `ALLOWED` au moment où le module est créé (étape 16).

### 21.2 La tension à nommer honnêtement

Firas demande que l'utilisateur n'ait pas à comprendre la complexité technique interne. C'est juste, et c'est faisable pour **la mécanique** : personne n'a besoin de savoir qu'un `PRES` est appliqué par `psfgen`, ni qu'un inventaire est relu dans le `.psf`.

Mais une interface ne peut pas rendre la **science** plus simple qu'elle n'est. Le champ de force, le modèle d'eau, le pH, l'objectif de simulation et le niveau de validation changent le résultat scientifique, et le principe P6 interdit d'en inventer une valeur par défaut. Masquer ces choix reviendrait à les prendre à la place de l'utilisateur, en silence.

La règle est donc : **la complexité d'implémentation est cachée, les décisions scientifiques ne le sont jamais.** Le rôle de l'interface est de les rendre **lisibles** — une par écran, formulée en français, avec sa conséquence —, pas de les supprimer. C'est aussi ce qui la rend utilisable par un non-programmeur : il n'a pas à écrire du YAML, il répond à des questions claires.

### 21.3 Technologie retenue : interface web servie localement

| Option | Pour | Contre |
|---|---|---|
| **Web servie en local** (`pepforge ui` ouvre le navigateur) | une seule base pour Windows et Linux, design moderne et responsive sans effort, s'installe par la même voie que le reste, suivi naturel des traitements longs (journal et progression en flux), déployable plus tard sur un serveur de laboratoire | le navigateur ne voit pas les chemins du disque : PepForge doit servir son propre explorateur de fichiers |
| **Bureau Qt** (PySide6) | boîtes de dialogue de fichiers natives, aucune dépendance au navigateur | empaquetage plus lourd sur deux systèmes, et l'aspect « moderne et épuré » demande nettement plus de travail |
| **Tkinter** | inclus dans Python | ne répond pas à l'exigence de design |

**Retenu : interface web servie localement.** L'exigence de responsive design que Firas mentionne la désigne, et c'est la seule option qui donne un résultat identique sur les deux systèmes sans double empaquetage. Le serveur ne quitte jamais la machine : il écoute sur l'interface locale uniquement.

### 21.4 Navigation : le pipeline est le modèle mental

L'interface suit l'ordre du traitement, parce que c'est l'ordre dans lequel l'utilisateur pense. Un assistant linéaire, avec les étapes visibles en permanence et la possibilité de revenir en arrière :

1. **Entrée** — séquence, fichier de structure, ou manifeste pour un lot (section 20).
2. **Peptide** — séquence relue, résidus affichés un par un, extrémités, pH, charge nette recalculée en direct.
3. **Substitutions aza-β³** — on clique une position, elle devient aza. C'est l'écran central de l'outil.
4. **Champ de force et moteur** — les choix, avec leur conséquence écrite.
5. **Solvatation** — boîte, sel, neutralisation.
6. **Vérification** — le niveau de validation 0 à 3, avec son coût annoncé avant de lancer.
7. **Résultats** — rapport, contrôles, fichiers.

À côté de l'assistant, trois espaces seulement : **Préparations** (l'historique, rejouable d'un clic), **Bibliothèque de résidus** (les entrées et leur provenance), **Outils** (les logiciels externes détectés, leur version, et ce qui manque).

### 21.5 Les écrans qui comptent

**Saisie de la séquence.** Le mode de lecture se choisit **avant** la saisie : code à une lettre ou code à trois lettres, deux boutons exclusifs dont aucun n'est choisi d'avance (section 5.7). Un bouton de bascule réécrit la saisie dans l'autre code (`convert_sequence`), comme le réglage `seq_view_format` de PyMOL, et change le mode en conséquence ; une conversion refusée laisse la saisie intacte et dit pourquoi. Le champ accepte la notation complète (`Ac-RGDS-NH2`, `Ac-Arg-Gly-Asp-Ser-NH2`). Validation en direct, caractère fautif surligné. Sous le champ, les résidus apparaissent en pastilles numérotées ; cliquer une pastille ouvre ce qu'on peut faire à cette position, substitution aza-β³ comprise. Les extrémités sont deux listes déroulantes, jamais un texte libre. La charge nette au pH choisi s'affiche en permanence et change en direct : c'est le retour qui apprend le plus vite à l'utilisateur.

**Les décisions ouvertes, visibles.** `PrepConfig.missing_decisions()` et `SolvationSpec.missing_decisions()` existent déjà et retournent la liste de ce qui reste à décider. L'interface les affiche comme une liste de contrôle permanente, et le bouton de lancement reste désactivé tant qu'elle n'est pas vide, avec le motif écrit. L'utilisateur ne peut donc pas lancer une préparation incomplète sans le savoir.

**Messages d'erreur.** Trois règles. Jamais de trace Python à l'écran. Chaque exception de `pepforge.exceptions` est traduite en trois éléments : ce qui s'est passé, où (le champ ou le fichier concerné), et **quoi faire**. Et pour un outil externe manquant, le message nomme l'outil, donne sa version attendue et renvoie à la procédure d'installation — pas « FileNotFoundError ».

**Résultats.** Le rapport de vérification en tableau, groupé par niveau, avec trois états distincts et visuellement distincts : réussi, échoué, **non vérifié**. Ce troisième état ne doit jamais ressembler au premier (section 18.9). Les fichiers produits sont listés avec leur rôle en français, pas seulement leur extension, et le YAML de configuration est téléchargeable pour rejouer la préparation.

**Mode lot.** Le manifeste est déposé, un tableau montre une ligne par peptide et son état, et une ligne en échec reste consultable sans interrompre les autres.

### 21.6 Identité visuelle

Sobre et professionnelle, au service de la lisibilité : fond clair, une seule couleur d'accent, le rouge et le vert réservés aux états de contrôle et à rien d'autre, polices système, espacement généreux, aucune ombre ni dégradé décoratif. Les tableaux et les nombres sont l'essentiel du contenu : ils passent avant tout effet. Un état « non vérifié » se distingue aussi par sa forme et son texte, et non par la seule couleur, pour rester lisible en cas de daltonisme et à l'impression.

**Progressive disclosure.** Deux niveaux : l'essentiel visible, l'avancé replié. Rien de ce qui change le résultat scientifique n'est replié, et tout ce qui est replié apparaît malgré tout dans le rapport final. C'est ainsi qu'on évite l'interface surchargée sans rien cacher.

### 21.7 Quand la construire

**Après que la chaîne fonctionne de bout en bout** (étape 10), sinon on habille un outil qui ne produit rien. L'interface est l'étape 16 du plan. En revanche, rien n'empêche de valider le **design** dès maintenant sur une maquette, ce qui est sans risque et évite de découvrir tardivement un désaccord sur l'ergonomie.
