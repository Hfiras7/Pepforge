# Bilan de l'étape C : liaisons covalentes entre résidus d'une structure importée

> Complément du module 3 ([PLAN-CONSTRUCTION.md](../PLAN-CONSTRUCTION.md), « Étape C »). Bilan établi le 28 novembre 2026. Implémentation approuvée par Firas le 28 novembre 2026 à 9 h 24 UTC, d'après la proposition architecturale remise le 27 novembre. Branche `claude/etape-c-liaisons-covalentes`, PR distincte empilée sur la PR #3 (module 3), qui n'est pas modifiée. Rien n'est fusionné.

**Statut final : TERMINÉ MAIS NON VALIDÉ.** Le code, les tests et la validation sur structures réelles sont faits. Revue de Firas du 28 novembre 2026 : architecture satisfaisante, tolérance de 0,45 Å validée **provisoirement**, vérifications demandées sur la tolérance, 6WQA et 4ZHL. Elles sont faites et rapportées en section 9 ; la vérification sur les coordonnées des entrées originales n'a pu être faite qu'en partie (section 9.4).

---

## 1. Nom et numéro

**Étape C (complément du module 3) : connectivité covalente.** Elle couvre la lecture des `CONECT`, l'évaluation des liaisons déclarées, la détection des liaisons possibles et les décisions explicites sur les composants.

Modules touchés : `io` (le gros du travail), `residues` (critère sourcé), `cli` (deux options). `core.structure` n'est touché que pour `Connection`, une classe ajoutée au module 3.

Non modifiés, conformément aux consignes : `Peptide`, `PrepConfig` (`input_structure` garde son nom, sans migration), la protonation (module 2), `checks` (module 4), les modules 1 et 2. Les étapes A et B ne sont pas commencées.

## 2. Objectif et responsabilités

L'objectif est de ne perdre aucun composant lié de façon covalente à un peptide importé, sans jamais prendre une proximité pour une preuve. Décisions de Firas du 28 novembre :

- lire les `CONECT` du PDB ;
- exploiter `LINK` et `struct_conn` ;
- examiner les composants non standards, y compris ceux des autres chaînes ;
- détecter les composants proches même sans déclaration ;
- ne pas confondre proximité et liaison ;
- rapporter chaque liaison avec ses composants, ses atomes, sa chaîne, sa distance, sa source et sa certitude ;
- ne supprimer aucun composant en silence : une situation ambiguë demande une décision explicite.

## 3. Fonctionnalités développées

### 3.1 Deux natures d'information, deux types

- **`Connection`** (liaison **déclarée**) ne vient que du fichier. Elle gagne cinq champs, tous avec une valeur par défaut :
  - les conformations des deux atomes (`first_alt_loc`, `second_alt_loc`) ;
  - les enregistrements qui la déclarent (`records`) ;
  - la distance écrite par le fichier (`reported_distance`) ;
  - l'indicateur de copie de symétrie (`symmetry_mate`, d'après `asu` de gemmi).

  Les trois derniers n'entrent pas dans l'égalité, pour qu'une même liaison lue en PDB et en mmCIF reste égale.
- **`BondAssessment`** (`io/connectivity.py`) est l'évaluation d'une liaison. Ses champs :
  - les deux atomes (`AtomRef` : résidu, chaîne, atome, élément, conformation) ;
  - la source : `declared`, `geometric` ou `undetermined` pour un contact métallique ;
  - le statut : `accepted`, `possible`, `incompatible`, `insufficient_data` ou `not_assessed` ;
  - la distance mesurée et la distance limite ;
  - le motif, la nature, les enregistrements d'origine et la distance écrite par le fichier ;
  - un **niveau de certitude** : `confirmed`, `declared_unverified`, `contradicted`, `possible`, `overlap` ou `undetermined`.

  Une évaluation géométrique n'est **jamais** convertie en `Connection`.

### 3.2 Lecture des `CONECT`

- Les `CONECT` sont résolus par numéro de série **dans le modèle lu**. Dans un ensemble RMN, le même numéro peut désigner un atome différent dans chaque modèle (cas réel de 1LCD.pdb).
- Ils sont fusionnés avec les `LINK`/`SSBOND` qui désignent les mêmes atomes. La liaison garde alors les deux sources (`records = ("LINK", "CONECT")`).
- Les `CONECT` internes à un résidu ne sont pas des liaisons entre résidus : ils sont comptés au journal.
- Un `CONECT` non attribuable (numéro absent, ou ambigu) va dans `Structure.metadata["conect_unresolved"]`. Il fait refuser l'extraction, pas la lecture.
- À l'écriture PDB, une liaison connue par `CONECT` seul est réécrite en `CONECT`, avec les numéros de série d'origine. Elle n'est jamais transformée en `LINK`, ce qui lui prêterait une nature non déclarée.
- En mmCIF, qui n'a pas de `CONECT`, elle part en `struct_conn` de type inconnu.

### 3.3 Critère sourcé : `residues/data/covalent_geometry.yaml`

- **Fenêtre de liaison** : `0,63 Å ≤ d ≤ r1 + r2 + 0,45 Å`, chargée par `load_covalent_criteria()` dans un `CovalentCriteria` validé.
- **Rayons** : rayons covalents de Cordero *et al.*, *Dalton Trans.* 2008, 2832-2838, pour 22 non-métaux et métalloïdes, avec la valeur sp3 pour le carbone. Ils ont été relus dans deux copies concordantes : ASE 3.30.0 et RDKit.
- **Tolérance de 0,45 Å et distance minimale de 0,63 Å** : convention de perception des liaisons d'Open Babel (O'Boyle *et al.*, *J. Cheminform.* 2011, 3:33). La tolérance est la même dans OEChem.
- **Élément inconnu** : un élément sans rayon sourcé lève `ReferenceDataError`, jamais une valeur supposée.
- **Métaux** : blocs s, d et f, plus les métaux du bloc p. Ils sont listés à part : c'est un classement de PepForge, sans valeur numérique. Un contact métallique est rapporté, jamais promu en liaison covalente, et l'ion reste un ion (section 11.4 de l'architecture).

### 3.4 Évaluation des liaisons déclarées (`assess_declared`)

Les cas sont examinés dans cet ordre :

1. Un métal est `not_assessed`.
2. Une nature hors covalent, disulfure ou inconnue est `not_assessed`.
3. Une liaison vers une copie de symétrie est `insufficient_data` : la distance dans l'unité asymétrique ne veut rien dire.
4. Un atome absent est `insufficient_data`.
5. Sinon, toutes les paires de conformations compatibles d'occupation non nulle sont mesurées :
   - dans la fenêtre : `accepted` ;
   - hors de la fenêtre : `incompatible`, ce qui fait refuser l'extraction ;
   - paires en désaccord entre elles : `insufficient_data`.

### 3.5 Détection des liaisons possibles (`find_contacts`)

- La recherche se fait par grille. Atomes lourds seulement, eaux exclues.
- Les liaisons peptidiques C(i)–N(i+1) et les liaisons déjà déclarées sont exclues.
- Champ de recherche : la chaîne, ses composants et **tous les résidus non polymères des autres chaînes**.
- Le rattachement se fait de proche en proche, par liaisons déclarées puis par contacts. Un ligand collé à un glycane lui-même lié est donc vu.
- Une paire d'atomes donne une seule entrée : `possible` dans la fenêtre, `overlap` sous 0,63 Å, `undetermined` avec un métal.
- Un contact entre deux résidus de la chaîne, ou avec un résidu polymère d'une autre chaîne, est refusé sans décision.

### 3.6 Composants et décisions explicites

- **`possibly_bound`** : nouvelle catégorie pour un composant relié par la **seule géométrie**. Il est refusé par défaut.
- **Site de rattachement** : chaque composant relié en reçoit un : `n_terminal`, `c_terminal`, `backbone_internal`, `side_chain`, `bridging` (plusieurs résidus) ou `component_only` (relié seulement à d'autres composants).
- **Décisions par chaîne**, dans `ExtractionDecisions`, donc aussi dans `ChainDecisions` :
  - `components: {étiquette: exclude_as_unbound | remove_covalent}`. La première décision s'applique à un composant `possibly_bound`, la seconde à un composant déclaré lié.
  - `contacts: {paire: not_bonded | disulfide}`. `disulfide` n'est admis que pour un contact SG–SG entre deux cystéines de la chaîne ; il ajoute un pont au `Peptide` par les liens existants, sans nouveau champ.
- **Décision sans objet** : une décision qui ne s'applique à rien (composant absent, catégorie inadaptée, contact inexistant) est refusée (`ConfigurationError`).
- **Composants écartés** : ils sont inscrits dans `Peptide.metadata["structure_components_set_aside"]` et au journal.
- **Ordre des refus** : les refus de l'étape C s'ajoutent après ceux du module 3, dont l'ordre est conservé.

### 3.7 Rapport, journal et ligne de commande

- `ExtractionReport.bonds` contient toutes les évaluations, et `uncertain_bonds()` liste celles qui demandent une décision.
- `to_dict()` ajoute `bonds` et `attachment`. Le journal reçoit les liaisons et le critère covalent, avec ses sources.
- `pepforge inspect --chain A` accepte `--component ÉTIQUETTE=DÉCISION` et `--contact PAIRE=DÉCISION`, deux options répétables. Le résumé affiche les liaisons et le critère.

## 4. Fichiers touchés

| Fichier | Nature |
|---|---|
| `src/pepforge/io/connectivity.py` | nouveau : types et évaluation des liaisons |
| `src/pepforge/residues/data/covalent_geometry.yaml` | nouveau : rayons, tolérance, métaux, sources |
| `src/pepforge/residues/structure_components.py` | `CovalentCriteria`, `load_covalent_criteria` |
| `src/pepforge/io/structure_formats.py` | provenance des liaisons, `CONECT` (lecture, résolution, écriture) |
| `src/pepforge/io/structure.py` | résolution des `CONECT` à la lecture et à la relecture d'écriture |
| `src/pepforge/io/peptide_extraction.py` | connectivité, `possibly_bound`, sites, décisions |
| `src/pepforge/io/__init__.py` | exports |
| `src/pepforge/core/structure.py` | `Connection` : cinq champs de provenance avec valeur par défaut |
| `src/pepforge/cli/main.py` | `--component`, `--contact`, résumé |
| `tests/test_connectivity.py` | nouveau, 51 tests |
| `tests/test_peptide_extraction.py`, `tests/data/structures/capped_peptide.{pdb,cif}` | géométries d'essai corrigées ; 4ZHL U et 6WQA A restent refusées, aucun contact admis comme pont (section 9) |
| `validation/etape_c_connectivity.py`, `validation/results/etape_c_connectivity.yaml` | nouveau : validation sur structures réelles |
| `validation/module03_structures.py`, `validation/results/module03_structures.yaml` | motifs de refus attendus mis à jour (section 7) |
| `docs/ARCHITECTURE.md`, `docs/PLAN-CONSTRUCTION.md`, `CLAUDE.md` | contrat de `io` et de `Connection`, étape C, commandes |

Aucune dépendance ajoutée : gemmi, déjà retenu au module 3, fournit `conect_map` et `asu`.

## 5. Tests et résultats

Résultats du 28 novembre 2026, Python 3.13, gemmi 0.7.5 :

- **`pytest` : 709 tests passent**. Le module 3 en avait 621 : il y en a 86 dans `test_connectivity.py` (dont 35 ajoutés après la revue, section 9.1) et 2 de plus dans `test_peptide_extraction.py`.
- **`ruff check .` et `ruff format --check .`** : propres.
- **`mypy --python-version 3.13`** : aucune erreur sur 42 fichiers. Lancé sans cette option, le poste local échoue sur les stubs de numpy, et non sur le code de PepForge : la CI le confirmera.
- **Tests synthétiques demandés par Firas** :
  - liaison intra-chaîne (contact SG–SG non déclaré, décisions `not_bonded` et `disulfide`, refus de `disulfide` hors SG–SG) ;
  - liaisons inter-chaînes (composant d'une autre chaîne rattaché, contact avec un polymère d'une autre chaîne refusé) ;
  - chélateur (DOTA fictif sur Lys, déclaré par `LINK` ou par `CONECT`, ou seulement proche) ;
  - ligand non lié qui le reste ; ion et contact métallique ;
  - composant pontant et rattachement de proche en proche ;
  - conformations alternatives, copie de symétrie, recouvrement, élément sans rayon ;
  - `CONECT` par modèle et non attribuables ;
  - relecture après écriture en PDB et en mmCIF.
- **Mutants faits à la main** : 16 sur 16 tués (11 avant la revue, 5 ajoutés en section 9.1), chacun défaisant une règle :
  - `CONECT` ignorés ;
  - autres chaînes non examinées ;
  - `asu` ignoré ;
  - `possibly_bound` écarté en silence ;
  - pas de rattachement de proche en proche, par liaison déclarée puis par contact ;
  - tolérance ignorée ;
  - recouvrement non distingué ;
  - métal traité comme covalent ;
  - contact non déclaré accepté sans décision ;
  - `disulfide` sans contrôle SG–SG.

## 6. Validation scientifique et références

Les fichiers sont les 14 structures du module 3 (copies de la suite de tests de Biopython 1.88) : 10 lisibles en un seul modèle, plus 1LCD.pdb lu au modèle 1. Le script est `validation/etape_c_connectivity.py`. Résultat : **129 contrôles, 0 échec** (107 avant la revue). Le module 3 garde **136 contrôles, 0 échec** ; il en avait 137, voir la section 7.

- **Pas de faux positif systématique** : toutes les liaisons possibles trouvées par la géométrie sont listées une à une et comparées à une table établie à la main (`EXPECTED_CONTACTS`). Il n'y en a aucune hors de 6WQA et 4ZHL.
- **Liaisons déclarées** : toutes celles qui touchent une chaîne peptidique sont confirmées par les coordonnées.
- **Composants non liés** : ils restent distingués des composants covalents.
  - Les ligands non liés restent `unbound_ligand` : ZMA de 6WQA, à 2,834 Å, est au-dessus de la fenêtre ; les quatre ligands de 4CUP aussi.
  - Le zinc reste un ion.
  - Les glycanes de 7DDO restent des modifications de chaîne latérale.
- **Contrôles positifs** : on retire une information d'un fichier réel et on vérifie que rien n'est perdu.
  - 7DDO sans `LINK` : les glycanes restent rattachés, par `CONECT`.
  - 7DDO sans `LINK` ni `CONECT` : ils deviennent `possibly_bound` et l'extraction est refusée. `exclude_as_unbound` les laisse hors du peptide, et le peptide le dit.
  - *Retirés après la revue de Firas* : les contrôles qui admettaient les contacts de 6WQA et de 4ZHL comme ponts. Ils sont remplacés par ceux de la section 9.

**Constats réels sur ces fichiers**, que le module 3 ne voyait pas :

- **6WQA, chaîne A** : quatre contacts SG–SG non déclarés, à 2,02–2,03 Å, entre les cystéines 71–159, 74–146, 77–166 et 259–262. La distance est celle d'un pont disulfure. La copie de Biopython ne contient que `_atom_site` : elle ne déclare aucune liaison, ni ces ponts ni d'autres (section 9.2).
- **4ZHL, chaîne U** : deux contacts SG–SG non déclarés, 50–111 à 1,64 Å et 136–201 à 1,51 Å. C'est **nettement plus court** que la somme des rayons covalents du soufre (2,10 Å). Ces contacts sont dans la fenêtre (au-dessus de 0,63 Å), donc rapportés comme possibles, mais une telle distance n'est pas celle d'un pont normal. La validation du fichier lui-même les classe en `_pdbx_validate_close_contact` (section 9.3). Refus maintenu, aucun pont admis.
- **1A8O.pdb** : 8 `CONECT` ne peuvent pas être attribués, parce que la copie de Biopython a renuméroté les atomes. L'extraction de la chaîne A est donc refusée tant que le fichier n'est pas corrigé ; la version mmCIF n'est pas touchée.

Références :

- Cordero B. *et al.*, « Covalent radii revisited », *Dalton Trans.* 2008, 2832-2838, DOI 10.1039/B801115J.
- O'Boyle N. M. *et al.*, « Open Babel: An open chemical toolbox », *J. Cheminform.* 2011, 3:33, DOI 10.1186/1758-2946-3-33.
- Format PDB 3.3 (`CONECT`, `LINK`, `SSBOND`) et dictionnaire mmCIF (`struct_conn`, `ptnr*_symmetry`).

## 7. Problèmes et limitations

1. **Tolérance de 0,45 Å** : c'est une convention logicielle (Open Babel, OEChem), pas une constante physique. Elle a été retenue parce qu'une détection trop large coûte une question, alors qu'une détection trop étroite coûte un composant perdu. **Elle est à confirmer par Firas.**
2. **Rayons de Cordero** : ils ont été relus dans deux copies concordantes (ASE, RDKit). L'article lui-même, en accès payant, n'a pas été consulté.
3. **Changements de comportement du module 3**, voulus mais à valider :
   - 4ZHL U, acceptée au module 3, est maintenant refusée sans décision.
   - 6WQA A et 1A8O.pdb A gagnent chacune un motif de refus.
   - Les motifs attendus de `validation/module03_structures.py` ont été mis à jour en conséquence, avec un commentaire qui renvoie ici.
   - Le contrôle de séquence de 4ZHL U ne peut plus se faire, puisque la chaîne est refusée. C'est pourquoi le module 3 compte 136 contrôles au lieu de 137. Il n'est pas remplacé : le refaire demanderait d'admettre les contacts, ce que Firas a exclu.
   - Les géométries d'essai synthétiques ont dû être corrigées : un O(i) était à 1,79 Å de N(i+1), des SG étaient placés hors de toute géométrie physique, et deux chaînes se chevauchaient. La détection les signalait à juste titre.
4. **Provenance après écriture en mmCIF** : une liaison connue par `CONECT` seul revient en `struct_conn` de type inconnu. Sa nature est conservée, mais sa source affichée change.
5. **Ce que l'étape ne fait pas** :
   - Les hydrogènes ne sont pas examinés.
   - Une liaison vers une copie de symétrie n'est pas mesurée (`insufficient_data`).
   - Les conformations alternatives sont comparées paire par paire compatible, sans choix d'une conformation.
   - Un contact avec un polymère d'une autre chaîne est refusé, jamais traité comme un composant.
   - La nature d'un contact métallique n'est pas déterminée. C'est une limite assumée : la coordination relève de la section 11.4.
6. **`mypy`** a été vérifié localement avec `--python-version 3.13` seulement (voir la section 5).
7. **Performances** : elles n'ont pas été mesurées au-delà des 14 fichiers ; le plus gros compte 18 550 atomes.

## 8. Statut final

**TERMINÉ MAIS NON VALIDÉ.**

Points attendus de Firas :

- la tolérance de 0,45 Å ;
- l'interprétation des contacts de 4ZHL et 6WQA ;
- les changements de comportement du module 3 ;
- la fusion de la PR, qui attend sa validation explicite.

Les étapes A et B ne commencent pas sans son instruction.

## 9. Revue de Firas du 28 novembre 2026 : vérifications demandées

### 9.1 Tolérance de 0,45 Å

Validée provisoirement comme convention algorithmique de détection. Elle reste dans `covalent_geometry.yaml` et nulle part ailleurs : un test relit le fichier, et un autre vérifie que les modules de détection n'ont aucune constante décimale.

- **Bornes** : la fenêtre est fermée des deux côtés. Pour sept paires d'éléments (N–C, C–C, O–C, S–S, S–C, P–O, Se–C) :
  - une distance égale à la limite est une liaison possible, et le flottant suivant ne l'est plus ;
  - à 0,63 Å exactement, la liaison est possible ; juste en dessous, c'est un recouvrement.
  - Une liaison déclarée suit la même fenêtre : confirmée à la limite, contredite 0,001 Å au-delà et sous 0,63 Å.
- **Sensibilité** (tests synthétiques, tolérance de 0,30 à 0,60 Å) :
  - la détection s'élargit de façon monotone ;
  - un contact placé à r1 + r2 + 0,45 Å n'est vu qu'à partir de 0,45 Å ;
  - le chélateur d'essai, à 1,45 Å de NZ, reste `possibly_bound` sur toute la plage.
- **Sensibilité sur les 14 fichiers réels**, de 0,30 à 0,80 Å (`sensitivity` dans les résultats de validation) :
  - le résultat est **identique de 0,30 à 0,70 Å** pour toutes les chaînes : 4 contacts dans 6WQA A, 2 dans 4ZHL U, aucun ailleurs ;
  - les premiers contacts qui ne sont pas des liaisons apparaissent à 0,80 Å (2XHE A ; 4ZHL U passe de 2 à 4). Ce sont des paires voisines de la chaîne, comme O(i)–N(i+1), à environ r1 + r2 + 0,75 Å ;
  - le choix entre 0,30 et 0,70 Å ne change donc rien sur ces fichiers, et 0,45 Å garde environ 0,3 Å de marge avant les premières fausses détections. Ce constat ne vaut que pour ces 14 fichiers.
- **Une liaison géométrique n'est jamais confirmée** :
  - un balayage de 0,05 à 3,00 Å sur les sept paires ne donne jamais que `possible` ou `overlap` ;
  - après une décision `disulfide`, le rapport garde le contact en `geometric / possible`, et la structure ne gagne aucune `Connection`. Le pont n'existe que dans le peptide, par décision nommée.
- **Mutants ajoutés**, tous tués : fenêtre ouverte côté géométrie, fenêtre ouverte côté déclaré, minimum ouvert, géométrie notée « confirmée », tolérance doublée.

### 9.2 6WQA, chaîne A : quatre contacts laissés possibles

Les données structurales, relevées dans la copie de Biopython, sont dans le tableau ci-dessous. Le seuil de la fenêtre S–S est de 2,55 Å.

| Paire | Distance SG–SG | Conformation | Occupation | Facteur B |
|---|---|---|---|---|
| CYS71–CYS159 | 2,028 Å | aucune | 1,00 / 1,00 | 45,2 / 46,0 |
| CYS74–CYS146 | 2,016 Å | aucune | 1,00 / 1,00 | 40,0 / 39,9 |
| CYS77–CYS166 | 2,022 Å | aucune | 1,00 / 1,00 | 32,1 / 31,9 |
| CYS259–CYS262 | 2,019 Å | aucune | 1,00 / 1,00 | 38,2 / 38,3 |

- Aucune des huit cystéines n'a de conformation alternative ni d'occupation partielle, sur aucun atome. Le fichier contient 98 atomes en conformation alternative, tous ailleurs.
- **Résultat géométrique, pas déclaration expérimentale** : la copie de Biopython ne contient **que** `_atom_site`. On n'y trouve ni `struct_conn`, ni méthode, ni résolution. L'absence de déclaration ne dit donc rien de l'entrée de la PDB.
- **Entrée originale** : l'en-tête RCSB (`files.rcsb.org/header/6WQA.cif`) a été lu par un outil de lecture web. Il donne un dépôt le 28 avril 2020, des révisions 1.0, 1.1 et 1.2 (18 novembre 2020, 18 octobre 2023, 13 novembre 2024), et ne contient pas de `struct_conn`. Mais il ne contient pas non plus de méthode ni de résolution, et sa lecture s'arrête au milieu d'une ligne. Il est donc probablement incomplet : **la question de savoir si l'entrée déclare ces ponts n'est pas tranchée**.
- Les contacts restent `possible`, l'extraction reste refusée sans décision, et aucun test ni contrôle ne les admet comme ponts.

### 9.3 4ZHL, chaîne U : refus maintenu

- Les deux paires : CYS50–CYS111 à 1,637 Å, et CYS136–CYS201 à 1,515 Å. Aucune conformation alternative ; occupation de 1,00 ; facteurs B de 49,0 / 49,1 et 45,6 / 49,6. Le fichier n'a aucun atome en conformation alternative. Méthode : diffraction X à 2,06 Å, R = 0,220, Rfree = 0,262.
- **La validation de la PDB incluse dans le fichier tranche** : `_pdbx_validate_close_contact` range précisément ces deux paires dans les contacts trop courts entre atomes **non liés**, à 1,51 et 1,64 Å. Les distances mesurées par PepForge concordent à 0,005 Å près, ce qu'un contrôle de la validation vérifie.
- Les trois autres ponts de la chaîne U (42–58, 168–182, 191–220) sont, eux, déclarés dans `struct_conn`.
- **Entrée originale** : l'en-tête RCSB actuel (révision 1.2 du 16 octobre 2024, la copie de Biopython étant la version 4.0001 de 2015) a été lu en entier. Il donne la même méthode, la même résolution, les mêmes quatre `struct_conn` (dont P:CYS1–P:CYS10) et les mêmes trois lignes `_pdbx_validate_close_contact`, avec les distances 1,51 et 1,64 Å. Les deux paires n'y sont donc pas déclarées comme ponts.
- Conséquence : extraction de la chaîne U refusée ; aucun test ne l'accepte. Les tests du module 3 qui portent sur autre chose (numérotation, décisions par chaîne) travaillent sur une copie d'essai privée des deux SG en cause (`_without_sg`), décrite comme telle : aucun pont n'y est admis.

### 9.4 Ce qui n'a pas pu être vérifié

- **Coordonnées des entrées originales** : le réseau de l'environnement cloud refuse les serveurs de la PDB (`files.rcsb.org`, wwPDB, PDBe, PDBj : erreur 403 du proxy, le 28 novembre 2026). Le serveur de coordonnées `models.rcsb.org` refuse les lecteurs automatiques (robots.txt). Seuls les en-têtes RCSB ont pu être lus, et seulement par un outil de lecture web dont le texte a été recopié ; la comparaison n'a donc pas été refaite par PepForge sur le fichier lui-même.
- Pour 4ZHL, la concordance des distances 1,51 et 1,64 Å entre l'en-tête actuel et la copie de 2015 suggère des coordonnées inchangées, sans le prouver atome par atome. Pour 6WQA, rien n'a pu être recoupé.
- Pour lever ce point, il faut l'un des deux : Firas fournit les deux fichiers mmCIF téléchargés depuis la PDB, ou l'accès réseau à `files.rcsb.org` est ajouté à l'environnement cloud.
