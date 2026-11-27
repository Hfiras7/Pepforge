# Bilan du module 1 — Modèle de données enrichi

> Étape 1 de [PLAN-CONSTRUCTION.md](../PLAN-CONSTRUCTION.md). Bilan établi le 27 novembre 2026, selon le modèle de la méthode de travail imposée. Code et tests : commit `1aa6b0a`.

**Statut final : TERMINÉ MAIS NON VALIDÉ.** Le code, les tests et la documentation sont livrés et la CI est verte sur Linux et Windows. La validation attend la revue de Firas prévue par le plan (section 8 ci-dessous).

---

## 1. Nom et numéro

**Module 1 — Modèle de données enrichi.** Modules de code touchés : `core` et `io.sequence`, comme prévu ; `checks`, `exceptions` et `cli`, pour la compatibilité (section 4). Le module `residues`, cité par le plan, n'a pas eu besoin de changer.

## 2. Objectif et responsabilités

Que `core` puisse représenter tout ce dont les étapes suivantes ont besoin (états de protonation, ponts disulfure, cyclisation, modifications de chaîne latérale, résidus non standards dans une séquence), **sans connaître aucun champ de force**, et sans qu'une déclaration puisse se perdre en silence dans une étape qui ne la traite pas encore.

## 3. Fonctionnalités développées

### 3.1 Un `Peptide` immuable

- Les champs ne se réaffectent pas, les conteneurs sont des tuples, et les métadonnées sont en lecture seule. Une liste passée au constructeur est figée, et un dictionnaire de métadonnées est recopié : le modifier ensuite ne change pas le peptide.
- Les métadonnées (texte saisi, extrémités supposées) sont **exclues de l'égalité** : `Ac-RGDS-NH2` et `Ac-Arg-Gly-Asp-Ser-NH2` donnent le même peptide. L'état de protonation, lui, fait partie de l'identité.
- Le peptide se copie (`copy.deepcopy`), se sérialise (`pickle`) et se convertit (`dataclasses.asdict`). Ce n'est pas un détail : la solution la plus courte, `types.MappingProxyType`, casse ces trois opérations (vérifié), ce qui aurait bloqué plus tard le traitement par lots.

### 3.2 Opérations, toutes renvoyant un nouveau peptide

| Opération | Effet | Refus (`PeptideModelError`) |
|---|---|---|
| `residue_at(position)` | résidu à une position, **numérotée à partir de 1** | position hors de 1 à n |
| `mutate(position, residue)` | remplace un résidu ; l'état de protonation de l'ancien n'est pas transmis | position qui porte un pont disulfure ou une modification |
| `with_n_terminus`, `with_c_terminus` | change une extrémité, qui devient déclarée (la marque « valeur par défaut » disparaît) | — |
| `with_protonation(position, état)` | fixe ou efface (`None`) l'état de protonation | clé mal formée |
| `with_link`, `without_link` | déclare ou retire une liaison | doublon, liaison absente |
| `with_modification`, `without_modification` | déclare ou retire une modification | doublon, modification absente |

**Règle de mutation.** Une mutation est refusée à une position engagée dans un pont disulfure ou portant une modification, parce qu'elle remplacerait la chaîne latérale dont dépend la déclaration, qui deviendrait fausse sans que rien ne le signale. Il faut retirer la déclaration d'abord, ce qui laisse une trace dans le journal de l'appelant. La cyclisation tête-queue passe par le squelette : elle n'empêche pas de muter le premier ou le dernier résidu.

### 3.3 État de protonation

`Residue.protonation` est facultatif. C'est une clé interne neutre (`his_delta`, `his_epsilon`, `his_protonated`) : lettres minuscules, chiffres et soulignés. Cette forme interdit d'office les noms de champ de force tels qu'on les écrit (`HID`, `HSD`...). Le modèle ne vérifie que la forme ; le sens de l'état revient à l'étape 2.

### 3.4 Liaisons et modifications déclarées

- `Link(kind, first, second)` avec `LinkKind.DISULFIDE` (pont disulfure) et `LinkKind.HEAD_TO_TAIL` (cyclisation tête-queue). Les positions sont remises dans l'ordre croissant, et un résidu ne peut pas être lié à lui-même.
- `Modification(key, position)`, en clé interne neutre.
- Invariants du peptide : positions dans les bornes, aucun doublon, une cyclisation relie exactement la position 1 et la dernière, et une position n'appartient qu'à un seul pont disulfure. Liaisons et modifications sont rangées dans un ordre canonique, si bien que l'ordre de déclaration ne change pas l'égalité.

### 3.5 Notation des non standards entre crochets

`Ac-RG[Aib]S-NH2` en notation à une lettre, `Ac-Arg-Gly-[Aib]-Ser-NH2` en notation à trois lettres (les crochets y sont facultatifs). La règle de lecture complète, en six points, est dans la docstring de [`io/sequence.py`](../../src/pepforge/io/sequence.py) : c'est le premier point à relire. Les crochets mal formés sont refusés avec un message qui situe l'erreur, par exemple :

```text
Erreur : Crochet ouvert (caractère 6) et non refermé avant le tiret (caractère 10)
dans 'Ac-RG[Aib-NH2' : un code entre crochets ne contient jamais de tiret
```

### 3.6 Contrôles et inventaire

- Nouveau contrôle `peptide.links_chemistry`. Un pont disulfure exige deux cystéines : sur un résidu standard sans thiol, c'est un échec ; sur un résidu non standard, c'est un échec « non vérifiable ». Une cyclisation exige une amine N-terminale et un acide C-terminal libres, et des résidus d'extrémité que la bibliothèque décrit.
- Nouveau contrôle `peptide.modifications_verifiable` : toute modification déclarée échoue comme « non vérifiable », puisqu'il n'existe pas encore de bibliothèque de modifications.
- **Inventaire attendu** : les ponts disulfure viennent désormais des liaisons déclarées, au lieu d'un `()` codé en dur qui aurait menti. Une cyclisation ou une modification, qui n'ont pas de champ dans l'inventaire, font **refuser** sa construction (`UnsupportedFeatureError`) au lieu d'être omises.

## 4. Fichiers créés ou modifiés

| Fichier | Nature |
|---|---|
| `src/pepforge/core/peptide.py` | réécrit : `Peptide` immuable, `Residue.protonation`, `Link`, `LinkKind`, `Modification`, opérations, constantes `N_TERMINUS_DEFAULTED` et `C_TERMINUS_DEFAULTED` |
| `src/pepforge/core/structure.py` | docstring : plus aucun nom de champ de force |
| `src/pepforge/core/__init__.py` | nouveaux exports |
| `src/pepforge/io/sequence.py` | notation entre crochets, règle de lecture documentée, résidus en tuple |
| `src/pepforge/checks/peptide_checks.py` | deux contrôles ajoutés, clés de métadonnées en constantes |
| `src/pepforge/checks/inventory.py` | ponts disulfure tirés des liaisons, refus des déclarations non représentables |
| `src/pepforge/exceptions.py` | `PeptideModelError`, `UnsupportedFeatureError` |
| `src/pepforge/cli/main.py` | mêmes constantes, sans changement de comportement |
| `tests/test_peptide.py` | créé, 64 tests |
| `tests/test_sequence.py`, `tests/test_checks.py`, `tests/test_inventory.py` | tests ajoutés (21, 8 et 4) ; les tests existants sont inchangés |
| `docs/ARCHITECTURE.md` | sections 4, 5.1, 5.5, 5.7, 7.2 et 17.1 |
| `docs/PLAN-CONSTRUCTION.md` | étape 1 (statut, critères), étape 7 (critère de refus) |

**Pourquoi `checks` a changé alors que le plan ne le citait pas.** Sans ces ajustements, une liaison ou une modification déclarée dans le modèle aurait passé les contrôles d'entrée et l'inventaire attendu sans être vue : le module aurait lui-même introduit une erreur silencieuse. C'est la règle « préserver la compatibilité entre les modules », pas une anticipation du module suivant.

## 5. Tests effectués et résultats

| Vérification | Résultat |
|---|---|
| `pytest` | **209 réussis** : 112 d'avant l'étape, inchangés, et 97 ajoutés |
| `ruff check`, `ruff format --check`, `mypy --strict` | propres |
| CI, [exécution 45989901326](https://github.com/Hfiras7/Pepforge/actions/runs/45989901326) | **5 tâches vertes** : Linux (Python 3.11, 3.12, 3.13) et Windows (3.11, 3.13) |
| Anciens fichiers de tests sur le nouveau code | les 49 tests de 7 fichiers, **dans leur version d'avant l'étape**, réussissent |
| Défauts injectés volontairement, un à la fois | **9 sur 9 détectés** par le test désigné |

**Critères d'acceptation du plan.**

1. *Aucun nom de résidu de champ de force dans `core`.* Le test cherche les six noms exigés par le plan et seize autres, AMBER (`ACE`, `CYM`, `ASH`, `GLH`, `LYN`) et CHARMM (`HSD`, `HSE`, `HSP`, `NTER`, `CTER`, `CT2`, `CT3`, `DISU`, `ASPP`, `GLUP`, `LSN`). La recherche est faite en mot entier et sans tenir compte de la casse. Pour prouver que ce test n'est pas vide de sens, je l'ai passé sur l'ancienne version de `core` : il y trouve `NHE` et `NME` dans `peptide.py`, `NHE` et `HIE` dans `structure.py`.
2. *Une mutation renvoie un nouveau `Peptide`, l'original est inchangé* : test dédié, plus les mêmes vérifications pour les extrémités et la protonation.
3. *Les notations existantes donnent le même résultat qu'avant* : prouvé par les anciens fichiers de tests passés sur le nouveau code, sans modification.

**Les neuf défauts injectés** : tri canonique des liaisons supprimé ; mutation permise sur une cystéine pontée ; métadonnées rendues modifiables ; marque d'extrémité supposée conservée après un choix explicite ; cyclisation acceptée entre deux positions quelconques ; ponts disulfure oubliés dans l'inventaire ; cyclisation inventoriée en silence ; pont disulfure accepté sans thiol ; crochets non vérifiés. Chacun fait échouer le test qui le vise, puis le fichier est restauré et la suite complète repasse.

## 6. Validation scientifique et références

Le plan prévoit « V : aucune » pour cette étape, faute de chimie quantitative. Ce qui a été vérifié :

- **Nomenclature** : codes à trois et à une lettre de IUPAC-IUB, *Nomenclature and Symbolism for Amino Acids and Peptides*, Eur. J. Biochem. 138, 9–37 (1984), déjà cités dans `standard_amino_acids.yaml`, inchangé.
- **Pont disulfure** : il relie les soufres des chaînes latérales de deux cystéines (cystine). D'où la règle : un résidu standard sans thiol est un échec, un non standard est non vérifiable tant que la bibliothèque ne décrit pas ses groupes.
- **Cyclisation tête-queue** : liaison amide entre l'amine α du premier résidu et le carboxyle α du dernier. Elle consomme ces deux groupes, et elle est donc incompatible avec une extrémité acétylée ou amidée.
- **Notation entre crochets** : convention de HELM pour les monomères désignés par plusieurs caractères. Zhang T., Li H., Xi H., Stanton R. V., Rotstein S. H., « HELM: A Hierarchical Notation Language for Complex Biomolecule Structure Representation », *J. Chem. Inf. Model.* 52(10), 2796–2806 (2012), DOI 10.1021/ci3001925. La référence a été vérifiée par sa notice Wikipédia, et la convention des crochets dans la documentation de l'analyseur HELM de Schrödinger. La page de spécification de la Pistoia Alliance s'est affichée vide : **je ne l'ai pas vérifiée sur la spécification officielle.**
- **Aucune donnée scientifique ajoutée** : ni charge, ni pKa, ni paramètre. L'acide α-aminoisobutyrique (`AIB`) n'apparaît que dans une bibliothèque construite par les tests.

## 7. Problèmes connus et limitations

1. **Protonation** : seule la forme de la clé est vérifiée. Qu'un état ait un sens pour son résidu (`his_delta` sur une alanine) ne sera contrôlé qu'à l'étape 2.
2. **Bibliothèque par défaut** : elle ne contient que les vingt standards. La notation entre crochets répond donc « résidu inconnu » pour tout non standard tant que l'étape 11A n'a pas ajouté les entrées : `pepforge parse "Ac-RG[Aib]S-NH2"` se termine en erreur, avec un message clair et sans trace Python.
3. **Liaisons et modifications** : elles se déclarent par l'API Python seulement, pas encore dans la notation des séquences (pas de `cyclo(...)`), ni dans la CLI, ni dans la configuration YAML.
4. **Cyclisation et modifications** : l'inventaire attendu les refuse. Un peptide qui en porte **ne peut donc pas encore être préparé**, et c'est voulu tant qu'aucune étape ne les traite. Les ponts disulfure, eux, sont inventoriés ; leur présence dans le `.psf` se vérifiera à l'étape 7.
5. **Modifications** : sans bibliothèque, toute modification déclarée échoue au contrôle d'entrée comme « non vérifiable ».
6. **Liste des résidus à thiol** : elle est fermée (cystéine seule) et écrite dans `checks`. Elle devra passer dans les données de la bibliothèque le jour où un thiol non standard (homocystéine, pénicillamine) y entrera.
7. **Peptide cyclique** : il garde ses extrémités déclarées (amine libre, acide libre), qui désignent les groupes consommés. Leur représentation dans les fichiers produits sera fixée avec la prise en charge de la cyclisation.
8. **Checklist 04** : elle n'a pas d'item pour les liaisons hors squelette. Les deux nouveaux contrôles n'ont donc pas de `checklist_item`.
9. **Provenance** : après une mutation, `metadata["input"]` reste le texte saisi à l'origine. L'historique des opérations revient au journal de l'appelant, comme le prévoit le plan. Aucune étape du pipeline n'appelle encore ces opérations, si bien que cette journalisation n'est pas exercée dans ce module.
10. **Affichage** : `one_letter_sequence` écrit `X` pour un résidu sans code à une lettre, comme avant. Il n'existe pas encore d'écriture d'une séquence en notation entre crochets.
11. **Livraison** : le module est sur la PR #1, et pas sur une PR dédiée comme le veut le plan (« une étape = une PR »), parce que je ne crée pas d'autre branche sans ton accord.

## 8. Ce qu'il faut pour valider

La revue prévue par le plan (définition de « terminé », point 6) :

1. **Relire la règle de lecture des séquences** (six points, docstring de `io/sequence.py`), en particulier la syntaxe `[CODE]` et le fait que `ALA` se lit A-L-A alors que `Ala` se lit alanine.
2. **Vérifier que les clés internes de protonation sont neutres** : lettres minuscules, du type `his_delta`, `his_epsilon`, `his_protonated`. Les noms définitifs seront fixés à l'étape 2.

Deux points à trancher, sans urgence :

- ajouter à la checklist 04 un item « liaisons hors squelette déclarées présentes dans la topologie » (limitation 8) ;
- confirmer que la règle de mutation (refus sur une position pontée ou modifiée, au lieu d'un retrait automatique de la déclaration) convient à ta façon de travailler.

Le statut passera à **TERMINÉ ET VALIDÉ** quand ces deux relectures seront faites.
