# Bilan du module 3 — Lecture et écriture de structures PDB / mmCIF, inventaire

> Étape 3 de [PLAN-CONSTRUCTION.md](../PLAN-CONSTRUCTION.md). Bilan établi le 27 novembre 2026, selon le modèle de la méthode de travail imposée. Livré sur une branche et une PR dédiées, partie de la branche par défaut (modules 1 et 2 fusionnés).

**Statut final : TERMINÉ MAIS NON VALIDÉ.** Le code, les tests, la validation sur 14 fichiers réels et la documentation sont livrés. Tous les contrôles passent. Le module n'est pas déclaré validé, pour trois raisons détaillées en section 8 : Firas n'a pas encore répondu sur la bibliothèque de lecture, les structures de validation ne sont pas celles de la thèse, et la revue prévue par le plan reste à faire.

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
| CI Linux (3.11, 3.12, 3.13) et Windows (3.11, 3.13) | voir la PR |

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

1. **Bibliothèque de lecture non confirmée par Firas.** Je lui ai posé la question au démarrage (carte de décision) et j'ai travaillé sur l'option recommandée, gemmi. Changer de bibliothèque ne toucherait que `structure_formats.py`.
2. **Les récepteurs de la thèse (DCC, TrkC, UNC5) ne sont pas dans la validation.** La PDB est inaccessible depuis l'environnement de développement. Le test V du plan prévoyait un récepteur ciblé : il reste à faire avec les fichiers ou liens de Firas (`pepforge inspect --url ...` sur son poste, ou un fichier local).
3. **Le test d'intégration sur une sortie de docking n'est pas fait** : aucune n'existe encore dans le projet.
4. **Lecture d'un seul modèle à la fois.** Un ensemble RMN se lit modèle par modèle ; l'analyse d'un ensemble entier n'est pas prévue avant l'étape 14.
5. **Numérotation de l'auteur seulement.** Les identifiants `label_asym_id` et `label_seq_id` du mmCIF ne sont pas gardés dans `Structure`. Ils servent à la relecture stricte, pas au-delà.
6. **Charge formelle non lue** (colonnes 79-80 du PDB, `pdbx_formal_charge`). Elle n'a pas encore d'usage, et la charge des ions sera traitée à la paramétrisation.
7. **Enregistrements `CONECT` et `ANISOU` non importés.** Les liaisons déclarées viennent de `SSBOND`, `LINK` et `struct_conn` ; `CONECT` (utile pour les ligands) n'est pas lu ni signalé.
8. **Le trou de numérotation n'est pas une preuve de rupture.** Il est signalé à l'inventaire et refusé par la correspondance avec un peptide (cas de la chymotrypsine 4ZHL). Le contrôle géométrique de continuité, qui tranchera, vient à l'étape 4.
9. **Une chaîne qui porte un résidu non observé est refusée par la correspondance**, même si ce résidu est à une extrémité (Aβ 1-16 dans 2BEG). C'est voulu (la séquence lue ne serait pas celle de la molécule), mais cela demandera une décision explicite (tronquer et coiffer) à l'étape 5.
10. **Un élément absent est déduit par gemmi du nom d'atome.** Le nombre d'atomes concernés est inscrit au journal, mais la déduction n'est pas refusée, parce que beaucoup de fichiers produits par des outils de modélisation n'ont pas de colonne d'élément.
11. **Le contenu entier est gardé en mémoire** (plafond par défaut 100 Mio), ce qui suffit pour un récepteur ou un complexe, pas pour un ribosome sans relever le plafond.
12. **Formats binaires non pris en charge** (BinaryCIF, MMTF) : refusés avec un message clair.

## 8. Ce qu'il faut pour valider

1. **Confirmer la bibliothèque** sur la carte de décision (gemmi + contrôles stricts recommandé).
2. **Relire les structures de validation choisies** (section 6.2) et les résultats (section 6.3).
3. **Relire les refus de la correspondance avec un peptide** (section 3.5 et limitations 8 et 9) : sont-ils ceux que tu veux, en particulier pour une chaîne dont les extrémités ne sont pas observées ?
4. **Revue prévue par le plan** : gestion des conformations alternatives et des codes d'insertion.
5. **Quand tu le souhaites**, une structure de tes récepteurs ou de tes peptides, par fichier ou par lien, pour refaire la section 6 sur un cas de la thèse.
