# Bilan du module 3 — Lecture et écriture de structures PDB / mmCIF, inventaire

> Étape 3 de [PLAN-CONSTRUCTION.md](../PLAN-CONSTRUCTION.md). Bilan établi le 27 novembre 2026, selon le modèle de la méthode de travail imposée. Livré sur une branche et une PR dédiées, partie de la branche par défaut (modules 1 et 2 fusionnés).

**Statut final : TERMINÉ ET VALIDÉ** (28 novembre 2026). Validation technique acceptée par Firas à l'audit final du 27 novembre ; fusion de la PR #3 autorisée par Firas le 28 novembre (11 h 07 UTC) et faite le même jour (commit `3e98201` de la branche par défaut). Historique des revues, conservé :

> Statut précédent : TERMINÉ MAIS NON VALIDÉ. Revue de Firas du 27 novembre 2026 : gemmi **validé**, structures de validation en **validation provisoire**, extraction d'un peptide **à corriger**. La correction est livrée (section 9) : rapport d'extraction par jonction, trois défauts corrigés, tests sur 4ZHL, 6WQA et 7DDO. L'extrémité non chimique décidée par Firas est livrée (section 9.6). Le module reste non validé tant que Firas n'a pas donné sa validation finale, requise avant la fusion de la PR #3.

---

## 1. Nom et numéro

**Module 3 — Lecture et écriture de structures PDB / mmCIF, et inventaire.** Modules de code touchés : `io` (le gros du travail) et `core.structure`, comme prévu par le plan. Trois autres modules sont touchés, avec justification en section 4 : `residues` (données sourcées sur les eaux et les groupes terminaux), `exceptions` (deux sous-classes de `StructureError`) et `cli` (commande `inspect`).

## 2. Objectif et responsabilités

Importer une structure 3D (expérimentale, modèle, ou plus tard sortie de docking) **sans la modifier**, depuis un fichier local ou une URL. En dresser l'inventaire et lister les décisions qu'elle demande. La réécrire sans perte, et établir le `Peptide` qu'elle décrit quand c'est possible. Rien n'est corrigé, rien n'est deviné : une structure qu'on ne peut pas lire exactement est refusée avec un motif.

Exigences ajoutées par Firas le 27 novembre 2026 :

- deux modes d'importation, fichier local et URL, avec un téléchargement contrôlé et une vérification du fichier obtenu ;
- une architecture générique, pour d'autres structures et d'autres formats ;
- à 16 h 25 UTC : ne choisir aucune structure à sa place ; puis à 16 h 35 UTC : choisir automatiquement des molécules de test pour valider le module.

## 3. Fonctionnalités développées

### 3.1 Obtention du fichier : `io/retrieval.py`

Les deux modes aboutissent au même objet, `RetrievedFile` : le contenu exact reçu, son origine, son nom, sa taille, l'empreinte **SHA-256 des octets reçus** (avant décompression, pour pouvoir la comparer à celle publiée par une source) et la date d'obtention.

- `load_local_file(chemin)` : fichier absent, dossier, fichier vide ou trop gros, refusés.
- `download_file(url, DownloadLimits(...))`, téléchargement contrôlé :

| Garde | Comportement |
|---|---|
| Schéma | `https` seulement ; `http` sur demande explicite (`allow_http=True`, `--allow-http`) ; `ftp`, `file` et les adresses sans hôte refusés |
| Redirections | suivies, mais aucune vers un autre schéma, et jamais de `https` vers `http` ; l'URL finale est gardée |
| Identifiants | une URL `https://nom:motdepasse@...` est refusée |
| Délai | `timeout_s` (60 s par défaut) pour chaque opération réseau |
| Taille | longueur annoncée vérifiée d'abord, puis lecture **arrêtée** au-delà de `max_bytes` (100 Mio par défaut), que l'annonce ait été honnête ou non |
| Réponse | code autre que 200 refusé (206 compris) ; page HTML refusée |
| Intégrité | longueur reçue comparée à la longueur annoncée (téléchargement tronqué) ; fichier vide refusé |
| gzip | décompressé sous le même plafond (protection contre une « bombe » gzip) ; flux tronqué ou suivi de données refusé |

Les plafonds sont des réglages techniques, pas des choix scientifiques : ils ont une valeur par défaut documentée, que l'appelant peut changer.

### 3.2 Lecture : `io/structure.py` et `io/structure_formats.py`

- `read_structure(chemin, journal=..., model=..., fmt=...)` et `read_structure_url(url, journal=..., model=..., fmt=..., limits=...)` passent par la même fonction, `parse_structure`.
- **Registre de formats** : `StructureFormat` (`pdb`, `mmcif`), les extensions reconnues, une signature de contenu, et deux registres `READERS` et `WRITERS`. Ajouter un format, c'est ajouter ses quatre éléments dans `structure_formats.py`, sans toucher au reste.
- **Format** : l'extension, le contenu et le format demandé doivent concorder. Un fichier sans extension (URL du type `.../entry?id=x`) est lu d'après son contenu, et le journal note que l'extension ne disait rien. Un contenu binaire (BinaryCIF, MMTF) ou non reconnu est refusé.
- **Bibliothèque** : gemmi 0.7.5 (MPL-2.0, roues Windows, annotations de type fournies), importé **uniquement** dans `structure_formats.py`, ce que vérifie un test d'architecture.
- **Lecture stricte**. J'ai constaté sur gemmi 0.7.5 que, sans rien signaler :
  - en PDB, la coordonnée `1x.104` est lue `1.0` ;
  - en PDB, une occupation absente devient 1,0 et un facteur B absent 20,0 ;
  - en mmCIF, `1x.0` et `?` deviennent NaN ;
  - un élément absent est déduit du nom de l'atome.
  
  PepForge relit donc lui-même chaque atome du fichier brut (colonnes de la spécification PDB 3.3, table `_atom_site` du mmCIF) avec des règles strictes, et en tire un **inventaire attendu**. La structure produite par gemmi est comparée à cet inventaire, atome par atome : pas un atome de plus ou de moins, même identité (modèle, chaîne, numéro, code d'insertion, résidu, nom, conformation), mêmes coordonnées. C'est la méthode de la section 17 de l'architecture, appliquée à la lecture. Une occupation ou un facteur B absents restent `None` : ils ne sont jamais remplis. Les atomes sans élément sont comptés et inscrits au journal.
- **Refus** : coordonnée, numéro ou occupation illisibles, ligne d'atome trop courte, atome en double, `MODEL` sans `ENDMDL`, fichier mmCIF à plusieurs blocs, fichier sans atome.
- **Modèles** : un fichier multi-modèles exige `model=` (`--model`) ; le premier n'est jamais pris en silence. Si l'en-tête `NUMMDL` annonce un autre nombre de modèles que le fichier n'en contient, le fichier est tronqué ou incohérent : refus.
- **Gardé tel quel** : conformations alternatives (avec occupations), codes d'insertion (`52` et `52A` sont deux résidus), hydrogènes, eaux, ions, ligands, indicateur `HETATM`, **liaisons déclarées** (`SSBOND`, `LINK`, `struct_conn` : disulfure, covalente, coordination métallique, liaison hydrogène) et **résidus déclarés non observés** (`REMARK 465`, `pdbx_unobs_or_zero_occ_residues` filtré sur les résidus de polymère non observés).
- **Journal** : une entrée `read_structure` avec la provenance (mode, chemin ou URL, URL finale, type de contenu, taille, SHA-256, date, gzip), le format, le format annoncé par l'extension, le modèle, le nombre de modèles, d'atomes, d'atomes sans élément ou sans occupation, de liaisons et de résidus non observés, et la version de gemmi.

### 3.3 Écriture : `write_structure(structure, chemin, fmt, journal=...)`

- PDB et mmCIF. Le fichier est écrit dans un fichier temporaire du même dossier, **relu depuis le disque**, comparé atome par atome à la structure (identité, coordonnées à l'arrondi du format près, élément, `HETATM`, occupation et facteur B), liaisons déclarées comprises. Il n'est placé à son chemin final qu'une fois la comparaison réussie ; sinon il est supprimé.
- Ce que le PDB ne peut pas porter lève `UnsupportedFeatureError` au lieu d'être tronqué : plus de 99 999 atomes, chaîne de deux caractères, numéro de résidu hors de [−999, 9 999], nom de résidu de plus de trois caractères, nom d'atome de plus de quatre, liaison hydrogène déclarée. Le message propose le mmCIF.
- Une occupation ou un facteur B inconnus (structure construite) sont écrits 1,00 et 0,00, comme le format l'exige, et leur nombre est inscrit au journal.

### 3.4 Inventaire : `io/structure_inventory.py`

`inventory(structure) -> StructureInventory`, sans modification. Par chaîne : atomes, acides aminés standards, résidus de type acide aminé, eaux, ions, composants non reconnus (par nom et nombre), premier et dernier résidu, trous de numérotation, codes d'insertion. Pour la structure : modèle lu et nombre de modèles, hydrogènes, résidus à conformations alternatives (avec la microhétérogénéité, quand le nom du résidu change selon la conformation), résidus déclarés non observés, liaisons déclarées par nature, atomes d'occupation nulle ou absente. Et la liste des **décisions à prendre** (`pending_decisions`) : modèle à confirmer, conformation à choisir, résidus à reconstruire, ruptures possibles, composants à garder ou retirer, ions, eaux, hydrogènes, liaisons à confirmer.

Règles de classement, toutes structurelles ou tirées d'un fichier sourcé :

- acide aminé standard : nom présent dans `standard_amino_acids.yaml` ;
- eau : `HOH` ou `DOD`, d'après `structure_components.yaml` (dictionnaire des composants chimiques de la wwPDB) ;
- ion : résidu d'un seul atome, qui n'est ni une eau, ni un acide aminé, ni un groupe terminal ;
- résidu de type acide aminé : standard, ou portant `N`, `CA` et `C` (ce qui compte `MSE` ou `AIB` dans la numérotation) ;
- tout le reste est un composant non reconnu, rapporté par son nom, jamais classé par supposition.

### 3.5 Correspondance : `structure_to_peptide(structure, chaîne, journal=...)`

> **Remplacé après la revue** : la version décrite ci-dessous est celle soumise à la revue. Elle confondait un saut de numérotation avec une rupture, et comptait ligands et glycanes comme résidus. Voir la section 9 pour l'extraction actuelle.

Établit le `Peptide` d'une chaîne, ou refuse avec le motif :

- résidu absent de la bibliothèque (`MSE`, `HYP`, nucléotide, ligand inclus dans la chaîne) : `UnknownResidueError` ;
- résidu déclaré non observé dans la chaîne, ou trou de numérotation : `StructureError`, parce que la séquence lue ne serait pas celle de la molécule ;
- microhétérogénéité : `StructureError` ;
- pont disulfure vers une autre chaîne, ou liaison covalente hors squelette : `UnsupportedFeatureError`. Une liaison peptidique déclarée (`LINK` autour d'un résidu `HETATM`) est reconnue et acceptée ; les liaisons hydrogène et de coordination ne font pas partie du peptide.

Les extrémités viennent de la structure quand elle les prouve : groupe `ACE`, `NH2` ou `NME`, atome `OXT` (carboxyle libre), hydrogènes `H1`/`H2`/`H3` (amine libre). Sans preuve, l'extrémité par défaut de la bibliothèque est appliquée **et marquée** (`N_TERMINUS_DEFAULTED`, `C_TERMINUS_DEFAULTED`), comme pour une séquence. Les ponts disulfure déclarés deviennent des `Link`. L'opération est journalisée.

### 3.6 Commande `pepforge inspect`

```bash
pepforge inspect structure.cif
pepforge inspect ensemble.pdb --model 1 --chain A
pepforge inspect --url https://exemple.org/structure.cif [--timeout 60] [--max-mb 100] [--allow-http]
```

Affiche en YAML la provenance, l'inventaire, les décisions à prendre et, avec `--chain`, le peptide. Une erreur s'affiche en une phrase, jamais en trace Python.

## 4. Fichiers créés ou modifiés

| Fichier | Nature |
|---|---|
| `src/pepforge/io/retrieval.py` | **nouveau** : fichier local et téléchargement contrôlé |
| `src/pepforge/io/structure_formats.py` | **nouveau** : registre de formats, lecture stricte, comparaison à gemmi, écriture ; seul fichier qui importe gemmi |
| `src/pepforge/io/structure.py` | **nouveau** : `read_structure`, `read_structure_url`, `parse_structure`, `write_structure` |
| `src/pepforge/io/structure_inventory.py` | **nouveau** : `inventory`, `structure_to_peptide` |
| `src/pepforge/io/__init__.py` | exports |
| `src/pepforge/core/structure.py` | `ResidueId`, `Connection`, `ConnectionKind` ; `Atom` reçoit code d'insertion, occupation, facteur B et refuse une coordonnée non finie ; `Structure` reçoit modèle, liaisons déclarées, résidus non observés, et des métadonnées en lecture seule. Tous les nouveaux champs ont une valeur par défaut : les constructions existantes restent valides |
| `src/pepforge/residues/data/structure_components.yaml` | **nouveau** : eaux, groupes terminaux, atomes qui prouvent une extrémité libre, avec source (CCD de la wwPDB, Westbrook 2015) |
| `src/pepforge/residues/structure_components.py` | **nouveau** : chargement et vérification du fichier |
| `src/pepforge/exceptions.py` | `StructureFormatError` et `RetrievalError`, sous-classes de `StructureError` |
| `src/pepforge/cli/main.py` | commande `inspect` |
| `pyproject.toml` | dépendance `gemmi>=0.7.5,<0.8` ; groupe `validation` (`biopython==1.88`) ; scripts de validation exemptés de la longueur de ligne |
| `.gitattributes` | fichiers d'essai gardés octet pour octet (pas de conversion des fins de ligne sous Windows) |
| `tests/data/structures/` | **nouveau** : `capped_peptide.pdb`, `capped_peptide.cif`, `nmr_two_models.pdb`, écrits à la main |
| `tests/test_structure_io.py`, `tests/test_retrieval.py`, `tests/test_structure_inventory.py` | **nouveaux** |
| `tests/test_cli.py`, `tests/test_architecture.py` | tests ajoutés (commande `inspect`, confinement de gemmi) |
| `validation/module03_structures.py`, `validation/module03_mutations.py` | **nouveaux** : validation sur structures réelles et défauts injectés |
| `validation/results/module03_structures.yaml` | **nouveau** : résultats complets de la validation |
| `docs/ARCHITECTURE.md` | sections 3 (bibliothèques tierces confinées), 4 (types), 5.7 (`io`) |
| `docs/PLAN-CONSTRUCTION.md` | étape 3 : statut, critères, test V corrigé |
| `README.md`, `CLAUDE.md` | commande `inspect`, ligne `io` |

Justification des modules touchés hors plan : `residues` parce que la règle « aucune donnée codée en dur » impose que les noms d'eau et de groupes terminaux viennent d'un fichier sourcé ; `exceptions` parce que chaque erreur a sa classe dédiée ; `cli` parce que Firas doit pouvoir importer ses structures sans écrire de Python. **La matrice de dépendances ne change pas** : `io` importait déjà `core`, `residues`, `report` et `exceptions`.

## 5. Tests effectués et résultats

Commandes : `ruff check .`, `ruff format --check .`, `mypy` (strict, sous Python 3.11) et `pytest`, toutes propres.

| Suite | Résultat |
|---|---|
| Tests existants (modules 1 et 2), fichiers inchangés | 460 sur 460 réussis sur le nouveau code |
| Nouveaux tests | 107 : 55 lecture et écriture, 26 téléchargement, 22 inventaire et correspondance, 3 commande `inspect`, 1 confinement de gemmi |
| **Total** | **567 réussis, 0 échec** |
| CI Linux (3.11, 3.12, 3.13) et Windows (3.11, 3.13) | **5 jobs sur 5 verts** sur le commit `56a8ed7` (PR #3, 27 novembre 2026), lint, typage et tests compris |

Tests sans Internet : fichiers écrits à la main, et un **serveur HTTP local** (127.0.0.1, port libre) qui simule chaque réponse à refuser (404, 206, HTML, trop gros annoncé, trop gros non annoncé, tronqué, trop lent, vide, redirections). La règle « jamais de `https` vers `http` » est vérifiée sur le gestionnaire de redirections, faute de serveur TLS de test.

**Défauts injectés** (`validation/module03_mutations.py`) : 33 gardes désactivées une à une (comparaison à gemmi, refus du double, plafonds, schéma, HTML, modèle non choisi, `NUMMDL`, format contredit, relecture de l'écriture, groupes terminaux, résidus non observés, trou de numérotation...). **33 sur 33 font échouer au moins un test.** La première passe en avait laissé 3 sans test (atome produit deux fois par le lecteur, dernier modèle sans `ENDMDL`, réponse 206) : les trois tests ont été ajoutés.

**Défaut trouvé pendant le développement** : le groupe amino terminal `NH2` n'a qu'un atome lourd. La première version de l'inventaire le comptait comme un ion, et la correspondance perdait l'extrémité amidée sans rien dire. Le test du peptide coiffé l'a montré ; les groupes terminaux sont désormais exclus des ions, et un défaut injecté vérifie que le test le détecte.

## 6. Validation scientifique et références

### 6.1 Ce qui a été comparé

Script [`validation/module03_structures.py`](../../validation/module03_structures.py), résultats complets dans [`validation/results/module03_structures.yaml`](../../validation/results/module03_structures.yaml). Chaque fichier est **importé par URL** (vrai téléchargement en `https`, avec SHA-256 au journal), puis comparé à :

1. **Biopython 1.88**, lecteur publié indépendant (Hamelryck et Manderick, *Bioinformatics* 19:2308, 2003, DOI 10.1093/bioinformatics/btg299), en mode strict pour le PDB : nombre de modèles, et pour **chaque modèle** l'ensemble des atomes (chaîne, numéro, code d'insertion, résidu, nom, conformation, coordonnées à 0,001 Å près), le nombre d'eaux et de résidus à conformations alternatives ;
2. **l'en-tête du fichier**, compté sans gemmi ni PepForge : résidus non observés et liaisons déclarées ;
3. **le même fichier dans l'autre format**, pour les trois entrées que la PDB fournit dans les deux ;
4. **l'écriture puis la relecture**, en PDB et en mmCIF ;
5. **la séquence déclarée** (`_entity_poly.pdbx_seq_one_letter_code_can`) pour les chaînes reconnues comme peptides.

### 6.2 D'où viennent les structures

Firas a demandé de choisir automatiquement les structures de test. Les serveurs de la PDB (`files.rcsb.org`, `data.rcsb.org`, `www.ebi.ac.uk`) sont **refusés par le réseau de l'environnement de développement** (HTTP 403 du proxy). J'ai donc pris des copies d'entrées de la PDB distribuées avec la suite de tests de Biopython 1.88, à l'étiquette Git `biopython-188` (adresse figée, `raw.githubusercontent.com`). J'ai vérifié que `1A8O.cif` de cette étiquette est identique, octet pour octet, à celui de l'archive source publiée sur PyPI (SHA-256 de l'archive conforme à celui affiché par PyPI). Ce ne sont **pas** des cibles de la thèse : ce sont des cas de lecture choisis pour couvrir chaque difficulté.

| Fichier | Ce qu'il exerce |
|---|---|
| 1A8O (PDB et mmCIF) | rayons X ; sélénométhionine `MSE` ; 88 eaux ; pont disulfure et `LINK` |
| 2XHE (PDB et mmCIF) | deux protéines, 6 315 atomes ; 142 résidus non observés |
| 1LCD (PDB et mmCIF) | RMN, 3 modèles ; protéine et ADN ; ion Na ; liaisons hydrogène déclarées en mmCIF |
| 2BEG (mmCIF) | RMN, 10 modèles ; fibrille Aβ(17-42), 5 chaînes ; 955 hydrogènes |
| 2BEG (PDB) | copie raccourcie par Biopython : 1 modèle, mais `NUMMDL 10` |
| 1AS5 (mmCIF) | RMN, 14 modèles ; peptide à hydroxyprolines `HYP`, extrémité `NH2` |
| 2n0n_M1 (PDB) | un modèle extrait d'un ensemble : `NUMMDL 20` |
| 7DDO (PDB) | cryo-ME ; conformations alternatives ; ion Zn ; glycanes `NAG` |
| 6WQA (mmCIF) | 98 atomes en conformations alternatives ; ligand |
| 4ZHL (mmCIF) | codes d'insertion (numérotation de la chymotrypsine) ; chaîne peptidique de 10 résidus |
| 4CUP (mmCIF) | conformations alternatives ; ligands ; résidus non observés |

### 6.3 Résultats

**112 contrôles, aucun écart.**

- **Atomes identiques à Biopython** dans les 12 fichiers lisibles, modèle par modèle : 56 892 atomes comparés au total, dont 18 550 sur les 10 modèles de 2BEG et 4 998 sur les 14 modèles de 1AS5.
- **Mêmes nombres** de modèles, d'eaux et de résidus à conformations alternatives que Biopython (1 pour 7DDO, 6 pour 6WQA, 2 pour 4CUP).
- **Résidus non observés et liaisons déclarées** égaux à l'en-tête : 142 pour 2XHE, 80 pour 2BEG, 15 pour 7DDO, 2 pour 4CUP ; 7 liaisons pour 1A8O, 14 pour 7DDO, 31 pour 1LCD en mmCIF.
- **PDB et mmCIF de la même entrée** : mêmes atomes pour 1A8O, 2XHE et 1LCD.
- **Écriture et relecture sans perte** en PDB et en mmCIF pour 23 des 24 cas. Le 24e est un **refus correct** : 1LCD en mmCIF déclare 27 liaisons hydrogène, que le PDB ne sait pas écrire. Le module refuse et propose le mmCIF au lieu de les perdre.
- **Fichiers tronqués refusés** : 2BEG.pdb et 2n0n_M1.pdb annoncent 10 et 20 modèles et n'en contiennent qu'un. **Biopython les lit sans erreur**, avec un seul modèle ; PepForge les refuse. C'est exactement l'erreur silencieuse que la section 17 interdit.
- **Correspondance avec un peptide** : la chaîne A de 1LCD (51 résidus) et la chaîne P de 4ZHL (`CPAYSRYIGC`, pont disulfure 1–10) donnent exactement la séquence déclarée par le fichier. Pour 1LCD, les deux extrémités sont prouvées par la structure (`H1` et `OXT`) ; pour 4ZHL, aucune ne l'est, et toutes deux sont marquées « par défaut ».
- **Refus motivés** : `MSE` (1A8O), `HYP` (1AS5), `NAG` (7DDO), `ZMA` (6WQA) et les nucléotides (1LCD) hors bibliothèque ; résidus non observés (2XHE, 2BEG, 7DDO, 4CUP) ; saut de numérotation 217 → 219 de la chymotrypsine (4ZHL, chaîne U).

### 6.4 Différences entre formats, constatées et non corrigées

Sur la même entrée, le PDB et le mmCIF ne portent pas tout à fait la même chose. Le module ne les « harmonise » pas, il les rapporte :

- 1A8O : `MSE` est en `HETATM` dans le PDB (120 enregistrements) et en `ATOM` dans le mmCIF (88). L'inventaire ne s'appuie donc pas sur l'indicateur `HETATM` pour classer un résidu ;
- 1LCD : le PDB déclare 4 liaisons (`LINK`), le mmCIF 31 (`struct_conn`, dont 27 liaisons hydrogène de l'ADN).

## 7. Problèmes connus et limitations

1. ~~Bibliothèque de lecture non confirmée par Firas.~~ *Validée le 27 novembre 2026 : gemmi, à condition que ses objets restent dans `io` (vérifié par `test_third_party_library_is_confined`).*
2. **Les récepteurs de la thèse (DCC, TrkC, UNC5) ne sont pas dans la validation.** La PDB est inaccessible depuis l'environnement de développement. Le test V du plan prévoyait un récepteur ciblé : il reste à faire avec les fichiers ou liens de Firas (`pepforge inspect --url ...` sur son poste, ou un fichier local).
3. **Le test d'intégration sur une sortie de docking n'est pas fait** : aucune n'existe encore dans le projet.
4. **Lecture d'un seul modèle à la fois.** Un ensemble RMN se lit modèle par modèle ; l'analyse d'un ensemble entier n'est pas prévue avant l'étape 14.
5. **Numérotation de l'auteur seulement.** Les identifiants `label_asym_id` et `label_seq_id` du mmCIF ne sont pas gardés dans `Structure`. Ils servent à la relecture stricte, pas au-delà.
6. **Charge formelle non lue** (colonnes 79-80 du PDB, `pdbx_formal_charge`). Elle n'a pas encore d'usage, et la charge des ions sera traitée à la paramétrisation.
7. **Enregistrements `CONECT` et `ANISOU` non importés.** Les liaisons déclarées viennent de `SSBOND`, `LINK` et `struct_conn` ; `CONECT` (utile pour les ligands) n'est pas lu ni signalé.
8. ~~Le trou de numérotation est refusé par la correspondance.~~ *Corrigé (section 9) : la liaison C–N est mesurée à chaque jonction ; un saut avec liaison est accepté et signalé.*
9. **Une chaîne dont des résidus ne sont pas observés est refusée par défaut**, à une extrémité comme dans une rupture (Aβ 1-16 dans 2BEG). Avec `observed_segment` ou `split`, chaque segment devient un peptide dont les bords portent l'extrémité non chimique `unobserved` (section 9.6) : à coiffer ou reconstruire explicitement avant la protonation.
10. **Un élément absent est déduit par gemmi du nom d'atome.** Le nombre d'atomes concernés est inscrit au journal, mais la déduction n'est pas refusée, parce que beaucoup de fichiers produits par des outils de modélisation n'ont pas de colonne d'élément.
11. **Le contenu entier est gardé en mémoire** (plafond par défaut 100 Mio), ce qui suffit pour un récepteur ou un complexe, pas pour un ribosome sans relever le plafond.
12. **Formats binaires non pris en charge** (BinaryCIF, MMTF) : refusés avec un message clair.

## 8. Ce qu'il faut pour valider

1. **Confirmer la bibliothèque** sur la carte de décision (gemmi + contrôles stricts recommandé).
2. **Relire les structures de validation choisies** (section 6.2) et les résultats (section 6.3).
3. **Relire les refus de la correspondance avec un peptide** (section 3.5 et limitations 8 et 9) : sont-ils ceux que tu veux, en particulier pour une chaîne dont les extrémités ne sont pas observées ?
4. **Revue prévue par le plan** : gestion des conformations alternatives et des codes d'insertion.
5. **Quand tu le souhaites**, une structure de tes récepteurs ou de tes peptides, par fichier ou par lien, pour refaire la section 6 sur un cas de la thèse.

## 9. Correction après la revue du 27 novembre 2026 : extraction d'un peptide

Orientations validées par Firas : gemmi encapsulé dans `io` ; seuil C–N de 1,8 Å, sourcé et modifiable, distance toujours rapportée ; classement explicite de la continuité, des sauts de numérotation, des ruptures et des composants non peptidiques ; ligands et glycanes jamais intégrés à la séquence ; mode strict par défaut, décisions `reject`, `split`, `observed_segment` journalisées ; aucune reconnexion ni reconstruction. Analyse préalable : [module-03-extraction-analyse.md](module-03-extraction-analyse.md).

### 9.1 Ce qui a changé

| Fichier | Changement |
|---|---|
| `src/pepforge/io/peptide_extraction.py` | **nouveau** : `analyze_chain`, `extract_peptides`, `structure_to_peptide` (forme stricte), `ExtractionReport`, `ExtractionDecisions` |
| `src/pepforge/io/structure_inventory.py` | l'ancienne correspondance en est retirée ; un saut de numérotation n'y est plus présenté comme une rupture |
| `src/pepforge/residues/data/structure_geometry.yaml` | **nouveau** : seuil 1,8 Å (Biopython 1.88, `PPBuilder`, Cock 2009, Hamelryck 2003), longueur de référence 1,329 Å (Engh et Huber 1991), validé par Firas |
| `src/pepforge/residues/structure_components.py` | `PeptideBondCriterion`, `load_peptide_bond_criterion(valeur)` ; une valeur imposée est inscrite comme telle et refusée si elle ne dépasse pas 1,329 Å |
| `src/pepforge/core/structure.py`, `io/structure_formats.py`, `io/structure.py` | `MissingAtom` ; atomes déclarés non observés lus (`REMARK 470`, `pdbx_unobs_or_zero_occ_atoms` avec `occupancy_flag` 1) et comptés au journal |
| `src/pepforge/exceptions.py` | `ExtractionRefusedError(StructureError)`, qui porte le rapport |
| `src/pepforge/cli/main.py` | `inspect --chain` imprime le rapport, même en cas de refus ; options `--breaks {reject,split}` et `--unobserved-termini {reject,observed_segment}` |
| `tests/test_peptide_extraction.py`, `tests/data/structures/reference/` | **nouveaux** : 39 tests ; 4ZHL, 6WQA, 7DDO compressés, empreinte SHA-256 vérifiée à chaque test |

La matrice de dépendances ne change pas.

### 9.2 Le rapport d'extraction

Chaque **jonction** porte les résidus, les atomes C et N utilisés (avec leur conformation, `A:SER6:C[A]`), chaque distance mesurée, l'écart de numérotation, la présence d'un code d'insertion, les résidus déclarés non observés entre les deux, la catégorie et la décision. Les conformations alternatives sont appariées par leur indicateur ; si elles ne concordent pas, la continuité est non mesurable.

| Catégorie de jonction | Critère | Par défaut |
|---|---|---|
| `continuous` | C–N ≤ 1,8 Å, numérotation n → n+1 sans code d'insertion | accepté |
| `numbering_discontinuity` | C–N ≤ 1,8 Å, saut, retour en arrière ou code d'insertion | **accepté et signalé** |
| `break_declared` | C–N > 1,8 Å, résidus déclarés entre les deux | refus, ou `split` |
| `break_undeclared` | C–N > 1,8 Å, aucun résidu déclaré | refus, ou `split` |
| `continuity_unmeasurable` | C ou N absent, conformations discordantes | refus |

| Catégorie de composant | Critère | Par défaut |
|---|---|---|
| `terminal_cap` | `ACE`, `NH2`, `NME`, à son extrémité | accepté (fixe l'extrémité) |
| `modified_residue` | N, CA, C présents, absent de la bibliothèque | refus, position donnée |
| `side_chain_modification` | relié à la chaîne par une liaison covalente déclarée, directement ou par un autre composant | refus, liaison donnée |
| `unbound_ligand` | aucune liaison covalente déclarée | exclu et listé, plus courte distance à la chaîne |
| `ion`, `water` | un atome ; nom d'eau | exclu et listé |

Aussi rapportés : extrémités non observées (`unobserved_termini`, refus ou `observed_segment`), résidus déclarés sans place dans la chaîne observée (contradiction, refus), atomes déclarés non observés, atomes de squelette absents (N, CA, C, O), segments.

### 9.3 Tests et résultats

- `ruff`, `ruff format`, `mypy --strict` propres ; **599 tests réussis, 0 échec** (39 nouveaux pour l'extraction, 3 pour les atomes déclarés manquants, 1 pour le seuil, 1 pour le rapport de refus en ligne de commande).
- Tests synthétiques à géométrie connue : chaque catégorie de jonction et de composant a un cas accepté et un cas refusé.
- Structures réelles (sans réseau) : 4ZHL chaîne U acceptée, 247 résidus, ARG217 → GLY219 à 1,316 Å classé `numbering_discontinuity`, 27 jonctions à code d'insertion, aucune distance au-delà de 1,8 Å ; 4ZHL chaîne P, 10 résidus, pont 1–10 ; 6WQA, `ZMA` classé ligand non lié (2,834 Å de la chaîne) et exclu, sauts +793 et −887 liés (1,355 et 1,351 Å), seule la rupture ALA1043 → ASP1060 (14,081 Å, non déclarée) bloque, et `split` donne deux segments de 254 et 137 résidus ; 7DDO, trois `NAG` classés modifications de chaîne latérale (Asn53, Asn90, Asn546) avec leur liaison `LINK`, Zn ion exclu avec ses trois coordinations, 596 jonctions continues ; 7DDO chaîne C, 14 résidus non observés à l'extrémité N et 1 à l'extrémité C.
- **Défauts injectés** (`validation/module03_mutations.py`) : 50 gardes désactivées une à une, dont 19 dans l'extraction, 1 sur le seuil et 2 dans la lecture des atomes déclarés manquants ; **50 sur 50 font échouer au moins un test**. La première passe sur l'extraction en avait laissé 2 sans test (appariement des conformations alternatives, liaison covalente vers une autre chaîne) : les deux tests ont été ajoutés.
- **Validation refaite** sur les 14 fichiers : **137 contrôles, 0 échec**. Les motifs de refus de chaque chaîne lisible sont comparés à ceux établis à la main dans l'analyse avant le code (`EXPECTED_REFUSALS`) : 20 chaînes conformes, et les 4 chaînes d'ADN de 1LCD refusées comme non peptidiques.

### 9.4 Défauts de la section 1 de l'analyse

| Défaut | Avant | Après |
|---|---|---|
| saut de numérotation pris pour une rupture (4ZHL U) | refus « chaîne interrompue » | extrait, saut signalé avec sa distance |
| ligand non lié compté comme résidu (6WQA `ZMA`) | refus « hors bibliothèque » | exclu et listé ; la vraie rupture est trouvée |
| glycane lié compté comme résidu (7DDO `NAG`) | refus « hors bibliothèque » | refus motivé : modification covalente de chaîne latérale, avec la liaison |
| rupture non déclarée invisible | non détectée | `break_undeclared`, par la géométrie |

### 9.5 Ce qui reste ouvert

1. ~~Extrémité non chimique~~ : décidée par Firas (option « nouvelle extrémité non observée ») et livrée, voir 9.6.
2. **Saut de numérotation avec liaison** (question A) et **granularité des décisions** (question B) : appliqués selon mes recommandations (accepté et signalé ; décision par catégorie et par chaîne), à confirmer.
3. Une liaison covalente non déclarée (fichier de docking sans `LINK`) ne peut pas être vue : le ligand est alors classé non lié. Sa plus courte distance à la chaîne est rapportée pour que l'utilisateur en juge ; aucun seuil n'est appliqué.
4. La liste complète des atomes absents par résidu, comparée au dictionnaire de la PDB, reste à l'étape 4.

### 9.6 Extrémité non chimique (décision de Firas, 27 novembre 2026, 19 h 33 UTC)

Option retenue : une **nouvelle extrémité « non observée »**, plutôt qu'une amine et un acide supposés.

- `core.peptide.Terminus` reçoit `chemical: bool = True` ; `termini.yaml` déclare `non_chemical_terminus` (`unobserved`, « non observée (la chaîne continue au-delà) »), sans alias et hors de `n_termini` / `c_termini` : une séquence ne peut pas la déclarer, et `ResidueLibrary.unobserved_terminus` la porte.
- `extract_peptides` la donne à tout bord de segment qui n'est pas une vraie extrémité (rupture coupée par `split`, résidus non observés au-delà avec `observed_segment`), même si le fichier y place un `OXT` ou un `H1`. Une vraie extrémité garde ses preuves (groupe terminal, `OXT`, `H1`), ou le défaut marqué. La nature des bords est gardée dans `Peptide.metadata["structure_edges"]` et au journal.
- **Module 2, modifié a minima** : `assign_protonation` et `theoretical_net_charge` refusent une extrémité non chimique (`UnsupportedFeatureError`, message qui dit de coiffer ou de reconstruire la chaîne d'abord). `format_sequence` refuse aussi de l'écrire. Tous les tests des modules 1 et 2 passent toujours.
- Tests ajoutés : 11 (extraction `observed_segment` et `split`, pont disulfure entre segments refusé, 6WQA découpé en 254 et 137 résidus, 7DDO C dont le glycane bloque toujours, refus de la protonation aux deux extrémités, refus d'écriture, chargement de la bibliothèque). **610 tests réussis, 0 échec.**
- Défauts injectés : 4 de plus sur les nouvelles gardes (refus de la protonation, refus d'écriture, bords N et C). **54 sur 54 détectés.** Une garde redondante a été retirée après qu'un défaut injecté a montré qu'elle ne changeait rien.
- Validation sur les 14 fichiers relancée : 137 contrôles, 0 échec (elle n'utilise que le mode strict).

### 9.7 Validation finale demandée par Firas (27 novembre 2026, 20 h 10 UTC)

Audit du code, des tests et du diff complet de la PR #3 (43 fichiers, aucun dans `checks/`, `build/`, `forcefields/`, `engines/` ni `pipeline/` : **module 4 non touché**). Contre la branche par défaut, aucune ligne de test n'est supprimée ; les seules lignes retirées hors du module 3 sont la signature interne `add_terminus` de la protonation, remplacée pour recevoir l'extrémité entière.

Noms des catégories : Firas a choisi de garder `numbering_discontinuity` et `unobserved` (carte du 27 novembre, 20 h 15 UTC), parce que la catégorie couvre aussi les codes d'insertion et les retours en arrière.

**Problèmes trouvés et corrigés**

| Problème | Conséquence | Correction |
|---|---|---|
| Le `Peptide` extrait ne gardait que le premier et le dernier identifiant du segment | la numérotation d'origine de chaque position (saut 217 → 219, codes d'insertion) n'était pas conservée dans le peptide | `Peptide.metadata["structure_residues"]` : identifiant d'origine de chaque position |
| Aucune représentation des décisions **par chaîne** | une configuration à plusieurs chaînes n'avait pas de forme, et rien n'empêchait d'appliquer à une chaîne la décision d'une autre | `ChainDecisions` (forme YAML, refus des chaînes, catégories et valeurs inconnues) et `extract_chains` |
| Message faux quand les atomes C et N existent mais dans des conformations différentes | le rapport annonçait « atome absent » avec une liste vide | message exact ; la jonction reste non mesurable et refusée |

**Tests ajoutés (11)** : numérotation d'origine de 4ZHL U conservée position par position ; preuve du saut (+2, 1,316 Å, atomes `U:ARG217:C` et `U:GLY219:N`) ; conformations sans paire compatible refusées ; deux chaînes à rupture dont seule A reçoit `split` (B reste refusée, journal des deux) ; deux chaînes avec chacune sa décision ; refus des chaînes, catégories et valeurs inconnues ; `split` seul ne lève pas le refus d'une extrémité non observée, et inversement ; les deux décisions ensemble, avec un `OXT` avant la rupture qui ne devient pas un acide libre ; segment extrait refusé par la protonation et le calcul de charge ; 4ZHL U et P extraites ensemble.

**Résultats**

- `ruff`, `ruff format`, `mypy --strict` propres ; **621 tests réussis, 0 échec**.
- Régression des modules 1 et 2 : les fichiers de tests de la branche par défaut (`test_peptide.py`, `test_sequence.py`, `test_library.py`, `test_protonation.py`, 287 tests) passent tels quels sur le nouveau code.
- Défauts injectés : **58 sur 58 détectés** (4 de plus : décision d'une chaîne propagée aux autres, chaîne non extraite acceptée, catégorie inconnue acceptée, message de conformations).
- Validation sur les 14 fichiers : **137 contrôles, 0 échec**, résultats identiques.

| Structure | Situation | Résultat obtenu |
|---|---|---|
| 4ZHL U | ARG217 → GLY219, C–N 1,316 Å | `numbering_discontinuity`, accepté et signalé ; 247 résidus, numérotation d'origine gardée |
| 6WQA A | ALA1043 → ASP1060, 14,081 Å | `break_undeclared` ; refus par défaut ; avec `split`, deux peptides de 254 et 137 résidus, bords `unobserved` |
| 6WQA A | `ZMA` | `unbound_ligand`, exclu et listé (2,834 Å) |
| 7DDO A | `NAG` sur Asn 53, 90, 546 | `side_chain_modification`, liaison `LINK` donnée, jamais un résidu ; refus motivé |

**Limites restantes**

1. Les atomes absents (hors C et N d'une jonction) sont rapportés sans décision : leur traitement est celui de l'étape 4.
2. Une liaison covalente non déclarée (fichier sans `LINK`) ne peut pas être vue ; le composant est classé non lié, avec sa plus courte distance à la chaîne.
3. Aucune borne basse n'est appliquée à la distance C–N : un chevauchement anormal (C–N très court) serait lu comme une liaison. Le contrôle des chevauchements est à l'étape 4.
4. Les décisions par chaîne existent dans l'API et sous forme YAML, mais aucune étape du pipeline ne les lit encore depuis `PrepConfig` ; la ligne de commande traite une chaîne à la fois.
5. Les récepteurs de la thèse (DCC, TrkC, UNC5) ne sont pas dans la validation, qui reste provisoire.
