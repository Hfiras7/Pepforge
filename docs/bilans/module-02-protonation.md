# Bilan du module 2 — Protonation et charge nette théorique

> Étape 2 de [PLAN-CONSTRUCTION.md](../PLAN-CONSTRUCTION.md). Bilan établi le 27 novembre 2026, selon le modèle de la méthode de travail imposée. Livré sur une branche et une PR dédiées, empilées sur la PR #1.

**Statut final : TERMINÉ MAIS NON VALIDÉ.** Le code, les tests et la documentation sont livrés. La source des pKa a été choisie par Firas. La validation attend trois choses de sa part (section 8) : relire les charges de référence, fournir le pKa de l'arginine lu sur l'article, et confirmer la fenêtre d'alerte. Tant que le pKa de l'arginine manque, **aucun peptide contenant une arginine ne passe cette étape**, ce qui inclut les peptides RGD.

---

## 1. Nom et numéro

**Module 2 — Protonation et charge nette théorique.** Modules de code touchés : `residues` (données et calcul), `checks` (contrôles de sortie), `pipeline` (étape `ProtonationStage`), comme prévu par le plan ; `core` et `config` pour deux ajouts justifiés en section 4 ; `exceptions` pour une nouvelle erreur.

## 2. Objectif et responsabilités

Choisir l'état de chaque groupe ionisable au pH de l'étude, calculer la charge nette qui servira de référence à tous les contrôles de charge ultérieurs (en particulier la charge relue dans le `.psf` à l'étape 7), et tracer chaque choix avec la source de son pKa. Sans jamais compléter une donnée manquante : un pKa absent, un résidu non décrit ou un tautomère indécidable arrêtent l'étape.

## 3. Fonctionnalités développées

### 3.1 Les données : `residues/data/ionizable_groups.yaml`

| Groupe | Porté par | pKa | États (charge) |
|---|---|---|---|
| amine α N-terminale | extrémité `free_amine` | 8,0 | `nterm_protonated` (+1), `nterm_deprotonated` (0) |
| carboxyle α C-terminal | extrémité `free_acid` | 3,7 | `cterm_protonated` (0), `cterm_deprotonated` (−1) |
| Asp | `ASP` | 3,9 | `asp_protonated` (0), `asp_deprotonated` (−1) |
| Glu | `GLU` | 4,3 | `glu_protonated` (0), `glu_deprotonated` (−1) |
| His | `HIS` | 6,5 | `his_protonated` (+1), `his_delta` (0), `his_epsilon` (0) |
| Cys (libre) | `CYS` | 8,6 | `cys_protonated` (0), `cys_deprotonated` (−1) |
| Tyr | `TYR` | 9,8 | `tyr_protonated` (0), `tyr_deprotonated` (−1) |
| Lys | `LYS` | 10,4 | `lys_protonated` (+1), `lys_deprotonated` (0) |
| Arg | `ARG` | **non renseigné** | `arg_protonated` (+1), `arg_deprotonated` (0) |

- **Noms définitifs des clés de protonation** (promis au module 1) : `<groupe>_protonated` et `<groupe>_deprotonated`, plus `his_delta` et `his_epsilon` pour les deux tautomères neutres de l'histidine (proton sur Nδ1 ou sur Nε2). Les trois clés d'histidine annoncées au module 1 sont conservées.
- Les 13 autres résidus standards sont classés **non ionisables** explicitement, comme les extrémités `acetyl`, `amide` et `n_methylamide`. Un test vérifie que les 20 résidus et les 5 extrémités de la bibliothèque sont tous classés : aucun n'est supposé neutre faute d'être cité.
- Le fichier est **validé au chargement** : chaque pKa cite une source déclarée, chaque groupe a une forme protonée et une forme déprotonée, les états d'une même forme ont la même charge, et la forme protonée porte exactement une charge de plus (un proton). 14 variantes fautives du fichier sont refusées par les tests.

### 3.2 Les calculs : `residues/protonation.py`

- `assign_protonation(peptide, ph) -> ProtonationResult`. Règle de choix, pour chaque groupe :
  1. un état déjà fixé (configuration ou peptide d'entrée) est **conservé et marqué « imposé »**, après vérification qu'il appartient au groupe ;
  2. sinon, la forme majoritaire au pH : protonée si pH < pKa, déprotonée si pH > pKa ;
  3. à pH = pKa exactement, **refus** : l'état doit être imposé ;
  4. si la forme majoritaire a deux états (histidine neutre), **refus** : le pKa ne départage pas les tautomères.
- `theoretical_net_charge(peptide, ph) -> float` : charge moyenne au pH, somme des charges moyennes de Henderson–Hasselbalch, la fraction protonée valant 1 / (1 + 10^(pH − pKa)).
- Deux charges sont donc produites : la charge **entière** des états retenus (`net_charge`), que la topologie devra porter, et la charge **moyenne** (`theoretical_net_charge`), qui dit à quel point le choix d'un état unique est une approximation. Pour H-GGG-OH à pH 7,4 : 0 et −0,2006.
- Une cystéine engagée dans un pont disulfure n'a plus de thiol : elle n'est pas ionisable.
- Refus explicites : pH hors de 0 à 14 (garde contre une faute de saisie), peptide cyclique et peptide modifié (`UnsupportedFeatureError`, la cyclisation consomme les extrémités, une modification peut ajouter un groupe ionisable), état fixé sur un résidu ou une extrémité sans groupe ionisable, ou sur une cystéine pontée (`ChemistryError`). Ce dernier point lève la **limitation 1 du module 1** : `his_delta` sur une alanine est désormais refusé.

### 3.3 L'étape : `pipeline.ProtonationStage`

- Lit le pH dans `PrepConfig.ph` ; s'il manque, `ConfigurationError`.
- Applique d'abord les états imposés de `PrepConfig.protonation` (site `"N"`, `"C"` ou position → état), puis choisit les autres.
- **Journal** : une entrée `impose_protonation` par état imposé, une entrée `assign_protonation` par groupe (site, résidu, groupe, état, charge, pH, pKa, **citation complète** de la source, raison du choix), et une entrée `net_charge`.

### 3.4 Les contrôles : `checks.PROTONATION_CHECKS`

| Contrôle | Issue | Ce qu'il attrape |
|---|---|---|
| `protonation.states_recorded` | `FAIL` | un état décidé absent du peptide, ou un état porté que rien n'a décidé |
| `protonation.net_charge` | `FAIL` | charge nette non entière, ou différente de la somme des états |
| `protonation.sources` | `FAIL` | un pKa sans source |
| `protonation.deduced_states_follow_ph` | `FAIL` | un état non imposé qui n'est pas la forme majoritaire |
| `protonation.near_pka` | `WARN` | \|pH − pKa\| < 1 : état à relire |
| `protonation.imposed_minority` | `WARN` | un état imposé qui n'est pas la forme majoritaire (choix légitime, mais à justifier) |

## 4. Fichiers créés ou modifiés

| Fichier | Nature |
|---|---|
| `src/pepforge/residues/data/ionizable_groups.yaml` | créé : groupes, états, pKa, sources, fenêtre d'alerte justifiée |
| `src/pepforge/residues/protonation.py` | créé : chargement validé, `assign_protonation`, `theoretical_net_charge` |
| `src/pepforge/core/protonation.py` | créé : `ProtonationResult`, `GroupAssignment`, `IonizationForm` (données seulement) |
| `src/pepforge/checks/protonation_checks.py` | créé : six contrôles |
| `src/pepforge/pipeline/stages.py` | `ProtonationStage` ajoutée |
| `src/pepforge/core/peptide.py` | `Terminus.protonation`, `with_n_terminus_protonation`, `with_c_terminus_protonation` |
| `src/pepforge/config.py` | champ `protonation` (états imposés), aller-retour YAML |
| `src/pepforge/exceptions.py` | `ReferenceDataError` |
| `src/pepforge/residues/library.py` | `_read_yaml` rendu public sous le nom `read_data_file` |
| `__init__.py` de `core`, `residues`, `checks`, `pipeline` | exports |
| `tests/test_protonation.py`, `tests/test_protonation_stage.py` | créés, 116 tests |
| `tests/reference/charges.yaml` | créé : 9 peptides de référence et 3 refus, avec le raisonnement |
| `tests/contracts.py`, `tests/test_interfaces.py` | contrat d'étape avec configuration facultative ; contrat de `ProtonationStage` |
| `docs/ARCHITECTURE.md` | sections 4, 5.1, 5.2, 5.4 et 6.4 |
| `docs/PLAN-CONSTRUCTION.md` | étape 2 : statut, prérequis, critères |

**Décisions d'architecture et leurs compromis.**

1. **`ProtonationResult` vit dans `core`, pas dans `residues`.** `checks` n'a pas le droit d'importer `residues` (matrice de dépendances). Mettre ce type de données dans `core` permet de le contrôler sans changer la matrice. Il ne contient aucune valeur : les pKa restent dans le fichier de données.
2. **`Terminus.protonation` ajouté au modèle du module 1.** L'état de l'amine N-terminale (chargée ou non) doit voyager avec le peptide jusqu'à l'adaptateur de champ de force, qui le traduira en patch. Le garder seulement dans le résultat aurait obligé chaque étape aval à transporter deux objets. Le changement est compatible : champ facultatif, `None` par défaut, et les 224 tests existants passent. Les deux nouvelles opérations **conservent la marque d'extrémité supposée**, contrairement à `with_n_terminus`, puisque l'extrémité elle-même ne change pas.
3. **`assign_protonation` renvoie un `ProtonationResult`, pas un `Peptide`** comme l'écrivait le plan : le peptide protoné est dans `result.peptide`, et le résultat porte en plus la trace de chaque choix, dont les contrôles ont besoin.
4. **`PrepConfig.protonation` ajouté**, vide par défaut, ce qui ne suppose rien. Sans lui, une histidine neutre ne pourrait pas passer l'étape depuis un fichier de configuration, et la règle 8 (relancer depuis la seule configuration) serait rompue.
5. **`build` n'intervient pas**, contrairement à ce qu'annonçait le tableau 6.4 de l'architecture : cette étape ne touche à aucune géométrie. Les hydrogènes de la structure 3D viendront avec l'étape 5.

## 5. Tests effectués et résultats

| Vérification | Résultat |
|---|---|
| `pytest` | **344 réussis** : 224 d'avant l'étape, et 120 ajoutés (116 dans les deux nouveaux fichiers, 1 contrat d'étape, 3 contrôles de dépendances pour les nouveaux fichiers source) |
| `ruff check`, `ruff format --check`, `mypy --strict` | propres |
| Anciens fichiers de tests sur le nouveau code | les fichiers de tests **dans leur version d'avant l'étape** réussissent sur le nouveau code (194 réussis, 1 ignoré hors du dépôt, le test d'architecture, qui lui est inchangé et passe dans le dépôt) |
| Défauts injectés volontairement, un à la fois | **14 sur 14 détectés** par le test désigné |
| CI GitHub (Linux 3.11 à 3.13, Windows 3.11 et 3.13) | sur la PR du module |

**Critères d'acceptation du plan.**

1. *Une extrémité acétylée ou amidée ne porte pas de charge ; une extrémité libre en porte une selon le pH* : `test_blocked_termini_carry_no_charge`, `test_n_methylamide_carries_no_charge`, et `test_free_termini_follow_the_ph` à pH 2, 7,4 et 11.
2. *Le pH vient de `PrepConfig.ph` ; s'il manque, `ConfigurationError`* : `test_missing_ph_raises_configuration_error`.
3. *Un pH proche d'un pKa produit un `WARN`* : `test_near_pka_is_a_warning_not_a_failure`, et la borne de la fenêtre est testée des deux côtés (écart 0,99 : alerte ; écart 1,00 : pas d'alerte).
4. *Le choix de chaque état est journalisé avec la source du pKa* : `test_each_choice_is_journaled_with_its_source`.

**Test d'intégration** : `SequenceInputStage` puis `ProtonationStage`, toutes deux par `run_stage` avec leurs contrôles (`test_sequence_then_protonation_through_run_stage`).

**Les quatorze défauts injectés** : signe de Henderson–Hasselbalch inversé ; pH = pKa tranché en silence ; tautomère d'histidine choisi en silence ; cystéine pontée traitée comme ionisable ; résidu non classé supposé neutre ; pKa manquant remplacé par une valeur ; état d'un autre groupe accepté ; contrôle « un proton entre les deux formes » affaibli ; cyclisation acceptée ; fenêtre d'alerte rendue inclusive ; contrôle des états portés toujours réussi ; états imposés de la configuration ignorés ; source du pKa absente du journal ; alerte de pKa transformée en succès. Chacun fait échouer le test qui le vise, puis le fichier est restauré et la suite complète repasse.

## 6. Validation scientifique et références

**Source des pKa, choisie par Firas le 27 novembre 2026** : les pKa mesurés sur des pentapeptides d'alanine par le groupe de Pace.

- Thurlkill R. L., Grimsley G. R., Scholtz J. M., Pace C. N., « pK values of the ionizable groups of proteins », *Protein Science* 15(5), 1214–1218 (2006), DOI 10.1110/ps.051840806.
- Grimsley G. R., Scholtz J. M., Pace C. N., « A summary of the measured pK values of the ionizable groups in folded proteins », *Protein Science* 18(1), 247–251 (2009), DOI 10.1002/pro.19.

**Ce que j'ai vérifié** : le tableau I de Grimsley 2009, lu sur PubMed Central (PMC2708032). Sa colonne « pK value in alanine pentapeptides » renvoie à Thurlkill 2006 et donne les huit valeurs du fichier. Sa note c précise que la valeur d'Asp (3,9) remplace celle de 3,67 publiée en 2006.

**Ce que je n'ai pas pu vérifier, et qui n'est donc pas dans le fichier** :

- **l'article de Thurlkill 2006 lui-même** : PubMed Central m'a renvoyé une page de vérification anti-robot, et le site de l'éditeur a refusé l'accès. Les conditions de mesure (séquences exactes, température, force ionique) restent à relever ;
- **le pKa de l'arginine** : absent du tableau de Grimsley 2009. La source prévue est Fitch C. A. *et al.*, « Arginine: Its pKa value revisited », *Protein Science* 24(5), 752–761 (2015), DOI 10.1002/pro.2647, dont je n'ai pu vérifier que l'existence et le DOI. Je ne reprends pas la valeur de mémoire.

**Pourquoi ce choix est défendable pour la thèse.** Ce sont des mesures faites sur des peptides courts, dans un même laboratoire, avec un même protocole : c'est le modèle le plus proche d'un peptide libre en solution. Pour le complexe peptide–récepteur, l'environnement déplace les pKa (écarts types de 0,5 à 2,7 unités selon le groupe dans le même tableau I) : ces valeurs y deviennent une première approximation, ce qui est écrit dans le fichier de données.

**Charges de référence calculées à la main** (`tests/reference/charges.yaml`), avec le raisonnement de chaque calcul écrit dans le fichier, indépendamment du code :

| Peptide | pH | Charge entière | Charge moyenne | À relire |
|---|---|---|---|---|
| H-GGG-OH | 7,4 | 0 | −0,2006 | extrémité N (écart 0,6) |
| Ac-GGG-NH2 | 7,4 | 0 | 0,0000 | — |
| Ac-KDEY-NH2 | 7,4 | −1 | −1,0039 | — |
| Ac-EEEE-NH2 | 7,4 | −4 | −3,9968 | — |
| H-GDG-OH | 2,0 | +1 | +0,9680 | — |
| Ac-GHG-NH2 | 5,0 | +1 | +0,9693 | — |
| Ac-GHG-NH2 (`his_epsilon` imposé) | 7,4 | 0 | +0,1118 | His (écart 0,9) |
| H-CG-OH | 7,4 | 0 | −0,2599 | extrémité N |
| H-CGC-OH, pont 1–3 | 7,4 | 0 | −0,2006 | extrémité N |

Refusés, comme attendu : Ac-RGDS-NH2 et H-RGDS-OH (pKa de l'arginine manquant ; charge entière attendue 0 dans les deux cas dès que ce pKa dépasse 8,4), et Ac-GHG-NH2 à pH 7,4 sans tautomère imposé.

**Comparaison à un logiciel de référence** : non faite. Les calculateurs de charge en ligne utilisent chacun leur propre jeu de pKa, si bien qu'un écart ne prouverait rien ; la vérification repose sur le calcul à la main avec les mêmes pKa. La comparaison pertinente viendra à l'étape 7 : la charge relue dans le `.psf` devra égaler `net_charge`.

## 7. Problèmes connus et limitations

1. **Arginine bloquée.** Tout peptide contenant une arginine est refusé tant que son pKa n'est pas saisi avec sa source. C'est voulu, mais cela bloque les peptides RGD de référence.
2. **pKa de groupes isolés.** Les interactions entre groupes voisins (quatre Glu adjacents, paire Asp–Lys) et l'environnement d'un complexe ne sont pas modélisés. Un calcul sur structure (PROPKA par exemple) pourra s'ajouter comme seconde source après l'étape 3, sans changer l'interface.
3. **Conditions de mesure non relevées** (Thurlkill 2006 non lu) : température et force ionique ne sont pas encore dans le fichier.
4. **Histidine neutre** : le tautomère doit toujours être imposé. C'est le bon comportement sans structure 3D, mais il demande à Firas une décision par histidine, à justifier dans `notes`.
5. **Fenêtre d'alerte** : un seuil de 1 unité, choix de prudence justifié par Henderson–Hasselbalch, à confirmer.
6. **Un état par groupe** : PepForge ne gère pas les états de protonation fractionnaires ni la dynamique à pH constant, hors du périmètre.
7. **Peptides cycliques et modifiés** : refusés par cette étape, comme par l'inventaire depuis le module 1.
8. **Résidus non standards** : aucun n'est encore classé. Chaque entrée de l'étape 11A devra déclarer ses groupes ionisables (ou leur absence) dans ce fichier ; un résidu aza-β³ ionisable aura sa propre entrée, puisque son environnement de squelette n'est pas celui d'un pentapeptide d'alanine.
9. **Interface** : la CLI (`pepforge parse`) n'affiche pas encore la charge ; ce n'était pas demandé par le plan pour cette étape.
10. **Inventaire attendu** : son champ `net_charge` reste à `None`. Il sera relié à `ProtonationResult.net_charge` à l'étape 7, qui en a besoin, et pas avant (pas d'anticipation).

## 8. Ce qu'il faut pour valider

1. **Relire les charges de référence** de la section 6 (et leur raisonnement dans `tests/reference/charges.yaml`).
2. **Fournir le pKa de l'arginine** lu sur Fitch *et al.* 2015 (valeur, incertitude, conditions), ou choisir une autre source publiée. Je l'ajouterai avec sa référence, et Ac-RGDS-NH2 passera dans les peptides de référence.
3. **Confirmer la fenêtre d'alerte** de 1 unité de pH.
4. Si possible, **relever les conditions de mesure** de Thurlkill 2006, pour les inscrire dans le fichier.

Le statut passera à **TERMINÉ ET VALIDÉ** quand les points 1 à 3 seront faits.
