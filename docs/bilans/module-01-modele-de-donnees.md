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
8. **Checklist 04** : elle n'avait pas d'item pour les liaisons hors squelette. *Résolu à la revue : item ajouté (section 9).*
9. **Provenance** : après une mutation, `metadata["input"]` reste le texte saisi à l'origine. L'historique des opérations revient au journal de l'appelant, comme le prévoit le plan. Aucune étape du pipeline n'appelle encore ces opérations, si bien que cette journalisation n'est pas exercée dans ce module.
10. **Affichage** : `one_letter_sequence` écrit `X` pour un résidu sans code à une lettre, comme avant. Il n'existe pas encore d'écriture d'une séquence en notation entre crochets. *Résolu au deuxième retour : `format_sequence` (section 10).*
11. **Livraison** : le module est sur la PR #1, et pas sur une PR dédiée comme le veut le plan (« une étape = une PR »), parce que je ne crée pas d'autre branche sans ton accord.

## 8. Ce qu'il faut pour valider

La revue prévue par le plan (définition de « terminé », point 6) :

1. **Relire la règle de lecture des séquences** (six points, docstring de `io/sequence.py`), en particulier la syntaxe `[CODE]` et le fait que `ALA` se lit A-L-A alors que `Ala` se lit alanine. *Remplacé au deuxième retour par le mode de lecture déclaré (section 10).*
2. **Vérifier que les clés internes de protonation sont neutres** : lettres minuscules, du type `his_delta`, `his_epsilon`, `his_protonated`. Les noms définitifs seront fixés à l'étape 2.

Deux points à trancher, sans urgence :

- ajouter à la checklist 04 un item « liaisons hors squelette déclarées présentes dans la topologie » (limitation 8) ;
- confirmer que la règle de mutation (refus sur une position pontée ou modifiée, au lieu d'un retrait automatique de la déclaration) convient à ta façon de travailler.

Le statut passera à **TERMINÉ ET VALIDÉ** quand ces deux relectures seront faites.

## 9. Revue de Firas, 27 novembre 2026

1. **Règle de lecture des séquences : approuvée.** Point soulevé à cette occasion : un bloc unique de trois lettres en majuscules qui est aussi un code (`ALA`) est lu en codes à une lettre. C'est le cas de **17 des 20 codes standards**, et `Ac-ALA-NHMe` (le dipeptide d'alanine, modèle classique en dynamique) se lisait **en silence** Ac-Ala-Leu-Ala-NHMe. La lecture ne change pas, mais elle est désormais **signalée** : `Peptide.metadata` (`AMBIGUOUS_READING`) dit ce qui a été lu et comment écrire l'autre lecture, et le contrôle `peptide.unambiguous_reading` produit un avertissement, affiché aussi par `pepforge parse`. Pour `GLU`, `LEU` et `PRO`, que la lecture en une lettre refuse, l'erreur indique d'écrire `Glu`, `Leu` ou `Pro`. *Remplacé au deuxième retour : le mode de lecture se déclare au départ, et cet avertissement est retiré (section 10).*
2. **Clés de protonation en minuscules : approuvées.**
3. **Item « liaisons hors squelette » dans la checklist 04 : ajouté**, sous l'identifiant `structure.links_and_modifications`, et relié à `peptide.links_chemistry`, `peptide.modifications_verifiable` et `inventory.disulfides`. Seule la copie du dépôt (`docs/04-checklist-qualite.md`) est modifiée ; le fichier d'origine du dossier d'instructions du projet ne l'est pas.
4. **Refus de muter une position pontée ou modifiée : conservé.** Un retrait automatique ferait disparaître une déclaration sans que l'utilisateur le voie, ce qui est exactement une erreur silencieuse. Le retrait reste une opération explicite (`without_link`, `without_modification`), que l'appelant journalise. Un balayage de positions (étape 11E) devra donc décider ouvertement quoi faire d'une cystéine pontée.

**Tests après la revue** : 224 réussis (15 ajoutés), et 3 défauts injectés sur 3 détectés (lecture ambiguë non signalée, erreur sans indication, avertissement transformé en succès).

**Statut : toujours TERMINÉ MAIS NON VALIDÉ**, en attendant que Firas confirme les ajouts des points 1 et 3, faits sur ma recommandation.

## 10. Deuxième retour de Firas, 27 novembre 2026 : mode de lecture déclaré et conversion

**Demande** : préciser dès le départ le mode de lecture, code à une lettre (monolitère) ou code à trois lettres (trilitère), et ajouter le passage d'une saisie d'un code à l'autre, comme dans PyMOL. Plus d'explications sur l'item ajouté à la checklist 04 (section 10.6).

### 10.1 Fonctionnalités

1. **Mode de lecture déclaré, jamais deviné.** `parse_sequence(texte, notation=...)` prend `SequenceNotation.ONE_LETTER` ou `SequenceNotation.THREE_LETTER` (ou leurs valeurs `one_letter`, `three_letter`). Le paramètre n'a **aucune valeur par défaut** : un appel sans mode est refusé. `ALA` est le tripeptide Ala-Leu-Ala en code à une lettre et l'alanine seule en code à trois lettres ; `Ac-ALA-NHMe` n'a plus de lecture devinée. Le mode est inscrit dans `Peptide.metadata` (`SEQUENCE_NOTATION`) et au journal de l'étape d'entrée.
2. **Règles propres à chaque code.** En code à une lettre, les résidus forment un seul bloc ; plusieurs blocs (`Arg-Gly-Asp`) sont refusés avec l'indication de choisir le code à trois lettres. En code à trois lettres, chaque résidu occupe un bloc, sans distinction de casse ; un bloc collé (`RGDS`) est refusé avec l'indication inverse. Quand un bloc refusé en code à une lettre est aussi un code à trois lettres (`Ala`, `GLU`), l'erreur le dit.
3. **Minuscules refusées en code à une lettre.** Dans l'usage, une minuscule désigne un acide aminé D (le `f` de `c(RGDfK)` est la D-phénylalanine). La bibliothèque ne décrit pas encore les acides aminés D : lire `f` comme une L-phénylalanine serait une erreur silencieuse, donc la lecture s'arrête et dit pourquoi. C'est un choix de ma part, au-delà de la demande.
4. **Conversion entre les deux codes.** `convert_sequence(texte, source, cible)` lit la séquence dans son code et la réécrit dans l'autre : `Ac-RGDS-NH2` ↔ `Ac-Arg-Gly-Asp-Ser-NH2`. `format_sequence(peptide, code)` écrit la forme canonique : résidus standards selon l'IUPAC (`R`, `Arg`), non standards entre crochets avec leur code de bibliothèque (`[AIB]`), extrémités sous leur forme usuelle, qui est le premier alias de `termini.yaml` (`H`, `Ac`, `OH`, `NH2`, `NHMe`). Une extrémité non déclarée **le reste** : convertir ne transforme pas une valeur supposée en valeur choisie.
5. **Relecture systématique.** Le texte produit est relu dans le code demandé et comparé au peptide d'origine, sur le principe de l'inventaire attendu contre sortie (architecture, section 17.1). Toute différence lève `UnsupportedFeatureError`. Une séquence écrite ne porte ni liaison hors squelette, ni modification, ni état de protonation : leur présence est refusée, pas omise.
6. **Ligne de commande.** `pepforge parse --notation 1|3` affiche le mode retenu et la séquence dans les deux codes. `pepforge convert --notation 1|3` réécrit la séquence dans l'autre code ; seul le résultat va sur la sortie standard, les avertissements sur la sortie d'erreur. Sans `--notation`, les deux commandes s'arrêtent sur un message qui dit quoi ajouter.
7. **Configuration.** `PrepConfig.sequence_notation` (`one_letter` ou `three_letter`) devient une décision exigée par `missing_decisions()` dès qu'une séquence est donnée, et toute autre valeur lève `ConfigurationError`. Sans ce champ, la configuration ne suffirait plus à relancer la même préparation (règle 8) : `ALA` aurait deux sens.
8. **Retrait de la lecture devinée.** `AMBIGUOUS_READING` et le contrôle `peptide.unambiguous_reading`, ajoutés à la première revue, sont supprimés : l'ambiguïté qu'ils signalaient ne peut plus se produire.

### 10.2 Fichiers touchés

- `src/pepforge/io/sequence.py` (réécrit : `SequenceNotation`, `SEQUENCE_NOTATION`, `parse_sequence`, `format_sequence`, `convert_sequence`), `src/pepforge/io/__init__.py`.
- `src/pepforge/residues/library.py` (forme écrite des extrémités : `n_symbols`, `c_symbols`, `n_terminus_symbol`, `c_terminus_symbol`) et `src/pepforge/residues/data/termini.yaml` (commentaire : le premier alias est la forme écrite ; aucune donnée changée).
- `src/pepforge/config.py` (`sequence_notation`, `SEQUENCE_NOTATIONS`). `config` ne peut importer que `exceptions` : la liste des modes y est recopiée, et un test vérifie qu'elle reste égale aux valeurs de `SequenceNotation`.
- `src/pepforge/pipeline/stages.py` (`SequenceInputStage(notation, ...)`, mode inscrit au journal), `src/pepforge/cli/main.py` (`--notation`, commande `convert`).
- `src/pepforge/checks/peptide_checks.py`, `src/pepforge/core/peptide.py`, `src/pepforge/core/__init__.py` : retrait de la lecture devinée.
- Tests : `test_sequence.py` (réécrit), `test_cli.py`, `test_config.py`, `test_library.py`, `test_checks.py` ; `test_peptide.py`, `test_inventory.py`, `test_interfaces.py` et `contracts.py` ne reçoivent que le paramètre `notation=`.
- Documentation : `docs/ARCHITECTURE.md` (P6, sections 5.7, 5.12, 20.2 et 21.5), `docs/PLAN-CONSTRUCTION.md` (étapes 1, 10, 11E bis et 16), `README.md`, `CLAUDE.md`.

### 10.3 Tests et résultats

- **333 tests réussis** (224 avant ce retour), `ruff` et `mypy --strict` propres.
- **Anciens tests** : leurs appels reçoivent seulement `notation=` (vérifié sur le diff, hors retours à la ligne de `ruff format`). Ont changé sur le fond, parce qu'ils encodaient la lecture devinée : `test_single_title_case_token_is_three_letter` et les trois tests de lecture ambiguë de `test_sequence.py`, `test_ambiguous_reading_warns_but_does_not_block` (`test_checks.py`), `test_parse_warns_on_ambiguous_reading` (`test_cli.py`), tous remplacés par des tests du mode déclaré. `test_complex_requires_receptor` reçoit `sequence_notation`, sans quoi le mode manquerait à la liste des décisions.
- **Aller-retour** : pour seize séquences, dont les vingt acides aminés standards à la suite, lire puis écrire dans chacun des deux codes redonne le même peptide.
- **Défauts injectés un à un : 16 sur 16 détectés**, chacun par le test qui le vise : minuscule lue comme un L ; plusieurs blocs acceptés en code à une lettre ; mode déclaré ignoré ; mode absent de `metadata` ; extrémité supposée écrite comme choisie ; relecture non comparée ; déclarations perdues sans refus ; code à trois lettres écrit en majuscules ; non standard écrit sans crochets ; indication vers l'autre code retirée ; extrémité écrite avec son dernier alias ; mode supposé par la CLI ; avertissements mêlés au résultat de `convert` ; mode non exigé par la configuration ; mode inconnu accepté par la configuration ; mode absent du journal. Caches de bytecode supprimés avant chaque essai.

### 10.4 Validation scientifique et références

- **PyMOL** : le réglage `seq_view_format` bascule l'affichage du visualiseur de séquence entre codes à une lettre (valeur 0, par défaut) et noms à trois lettres (valeur 1). Lu sur le miroir du PyMOL Wiki, <https://wiki.pymol.org/index.php?title=Seq_view_format>, page modifiée pour la dernière fois en 2014. PyMOL change un **affichage** ; PepForge réécrit le **texte saisi**. L'idée est la même, l'objet diffère.
- **Symboles des résidus** : IUPAC-IUB 1984, déjà cité dans `standard_amino_acids.yaml` ; l'écriture `Arg`, `Gly` en est la forme à trois lettres.
- **Minuscule pour un acide aminé D** : convention d'usage, illustrée par `c(RGDfK)`. **Je ne l'ai pas vérifiée sur un texte normatif** ; le refus ne dépend pas de ce point, puisqu'il suffit que la lecture soit douteuse pour être refusée.
- Aucune donnée scientifique ajoutée : ni charge, ni pKa, ni paramètre.

### 10.5 Problèmes et limitations

1. **Corps à une lettre écrit comme un alias d'extrémité N.** `AC-NH2` (Ala-Cys amidé, extrémité N non déclarée) se lit avec `AC` comme l'extrémité acétylée, et la lecture **échoue** sur `NH2` ; `H-NH2` (histidine amidée) de même. L'échec est bruyant, jamais silencieux, et déclarer l'extrémité N suffit (`H-AC-NH2`). Pour la même raison, `format_sequence` refuse d'écrire l'histidine amidée sans extrémité N en code à une lettre (`H-NH2` ne se relirait pas) ; en code à trois lettres, `His-NH2` passe.
2. **Acides aminés D** : non représentés. En code à une lettre, la minuscule est refusée ; en code à trois lettres, `D-Phe` échoue sur le bloc `D`, résidu inconnu.
3. **Forme canonique** : la conversion n'emploie pas l'orthographe saisie des extrémités. `NH2-RGDS-CO-NH2` devient `H-Arg-Gly-Asp-Ser-NH2`, même molécule, autre écriture.
4. **`PrepConfig.sequence_notation`** est exigé et sérialisé, mais le pipeline ne lit pas encore la configuration ; ce sera l'étape 10, où un critère d'acceptation est ajouté. La bascule de l'interface est décrite (architecture, section 21.5 ; plan, étape 16), pas construite.
5. **La CLI** accepte `1`, `3`, `one_letter` et `three_letter` ; l'API et la configuration, seulement `one_letter` et `three_letter`.

### 10.6 L'item « liaisons hors squelette » de la checklist 04, expliqué

**Ce qu'est la checklist 04.** `docs/04-checklist-qualite.md` est la liste des contrôles à passer avant toute simulation, avec une règle : un point non vérifié est un point en échec. Chaque contrôle automatique de PepForge y renvoie par son `checklist_item`. Le rapport peut ainsi dire, ligne par ligne, ce qui a été vérifié et par quel contrôle.

**Le trou.** Le module 1 a donné au modèle des déclarations qu'il n'avait pas : ponts disulfure, cyclisation tête-queue, modifications de chaîne latérale. Aucune ligne de la checklist ne les couvrait. Or le catalogue des défaillances silencieuses (architecture, section 17.2) cite justement « pont disulfure non déclaré » : si le patch du pont n'est pas appliqué, `psfgen` construit deux cystéines libres au lieu d'un pont, sans erreur et avec un code de retour nul.

**L'item ajouté**, dans la partie « Structure » : « Liaisons hors squelette (ponts disulfure, cyclisation) et modifications de chaîne latérale : chimiquement possibles, vérifiables, et présentes dans la topologie produite. » Trois contrôles y sont rattachés, à deux moments du pipeline :

- à l'entrée, `peptide.links_chemistry` : un pont disulfure ne relie que deux cystéines, une cyclisation tête-queue exige une amine libre et un acide libre aux extrémités ;
- à l'entrée, `peptide.modifications_verifiable` : une modification déclarée échoue tant qu'aucune bibliothèque de modifications ne permet de la vérifier ;
- à la sortie de la topologie (étape 7), `inventory.disulfides` : les ponts attendus, tirés du `Peptide`, sont comparés aux liaisons S–S réellement présentes dans le `.psf`.

**Exemple.** `Ac-CRGDC-NH2` avec un pont 1–5. Si le patch est oublié, le `.psf` contient zéro pont au lieu d'un : `inventory.disulfides` échoue, et le rapport range cet échec sous cette ligne de la checklist. Sans l'item, l'échec aurait bien été détecté, mais il n'aurait correspondu à aucune ligne : la checklist, qui est ce qu'on relit avant de lancer une simulation, aurait paru complète alors qu'un de ses aspects n'y figurait pas.

**Portée de la modification** : seule la copie du dépôt (`docs/04-checklist-qualite.md`) est modifiée, pas le fichier d'origine du dossier d'instructions du projet.

### 10.7 Statut

**TERMINÉ MAIS NON VALIDÉ.** Le statut passera à **TERMINÉ ET VALIDÉ** quand Firas aura confirmé ces ajouts (mode de lecture, conversion, refus des minuscules) et l'item de la checklist 04.
