# PepForge — Plan de construction

> Construction de PepForge en étapes courtes, chacune livrée par une PR séparée, revue et testée avant de passer à la suivante. L'architecture cible est décrite dans [ARCHITECTURE.md](ARCHITECTURE.md).

## Mode d'emploi

- **Une étape = une PR** (ou quelques PR si l'étape est longue), sur une branche dédiée, fusionnée seulement quand tous ses critères d'acceptation sont cochés.
- **Définition de « terminé »**, commune à toutes les étapes :
  1. CI verte : `ruff check`, `ruff format --check`, `mypy --strict`, `pytest` sur Python 3.11, 3.12 et 3.13 ;
  2. `tests/test_architecture.py` passe : aucune dépendance interdite ;
  3. toute nouvelle donnée scientifique a sa source (publication, base de données, version) ;
  4. toute opération qui modifie une structure est journalisée ;
  5. `ARCHITECTURE.md` est à jour si un contrat ou un type d'échange a changé ;
  6. le modèle de PR est rempli, et Firas a relu la partie scientifique.
- **Étapes bloquées par une décision** : elles sont marquées 🔒. On peut avancer sur les étapes suivantes qui n'en dépendent pas.
- **Tests** : quatre niveaux (voir [ARCHITECTURE.md, section 7](ARCHITECTURE.md#7-vérification--revue-contrôles-et-tests)). Abréviations ci-dessous : **U** unitaire, **K** contrat, **I** intégration, **V** validation scientifique.

> **Fil central, ajouté le 26 novembre 2026.** La fonction centrale de PepForge est la substitution d'un acide aminé par un résidu **aza-β³**, pour n'importe quel acide aminé et n'importe quelle position, suivie de la génération des fichiers de dynamique moléculaire (voir [ARCHITECTURE.md, section 10](ARCHITECTURE.md#10-cas-dusage-central--substituer-un-acide-aminé-par-un-résidu-aza-β³)). Les étapes 11A à 11E ci-dessous portent ce fil. Elles ne sont plus une extension tardive : dès que l'étape 4 est passée, le fil aza-β³ avance **en parallèle** des étapes 5 à 10, parce que la dérivation des paramètres est l'étape la plus longue et la plus risquée du projet. L'étape 10 (bout en bout) doit être démontrée sur un **pseudopeptide aza-β³**, et pas seulement sur un peptide standard.

```mermaid
flowchart LR
    E0[0. Socle] --> E1[1. Modèle enrichi]
    E1 --> E2[2. Protonation et charge]
    E1 --> E3[3. PDB / mmCIF]
    E3 --> E4[4. Contrôles structure]
    E2 & E4 --> E5["5. Construction 3D 🔒"]
    E2 --> E6[6. Adaptateur CHARMM36m]
    E5 & E6 --> E7["7. Topologie NAMD"]
    E7 --> E8[8. Solvatation]
    E8 --> E9[9. Protocole]
    E9 --> E10[10. Rapport et peptide libre de bout en bout]
    E4 --> E11["Fil aza-β³ : 11A à 11E<br/>fonction centrale"]
    E11 --> E10
    E10 --> E12[12. Complexe peptide–récepteur]
    E10 --> E13[13. Second moteur et conversion]
    E10 --> E14[14. Analyses]
```

---

## Étape 0 — Socle d'architecture ✅ (cette PR)

**Objectif** : fixer les modules, leurs contrats et les règles de vérification avant d'écrire les fonctionnalités.

**Livrables** : `docs/ARCHITECTURE.md`, ce plan, `core.Structure`, module `checks`, module `build` (interfaces), `EngineWriter` revu (`write_topology`, `solvate`, `EngineFiles`), `ParameterSet` enrichi, `pipeline.Stage` / `run_stage`, `ExternalTool` déplacé dans `tools`, `SolvationSpec`, `tests/test_architecture.py`, `tests/contracts.py`, modèle de PR.

**Critères d'acceptation**

- [x] Chaque module a un rôle, des entrées, des sorties et des interdits écrits.
- [x] La matrice de dépendances est vérifiée automatiquement.
- [x] Chaque interface a une fonction de contrat réutilisable.
- [x] Une étape réelle (`SequenceInputStage`) passe par `run_stage` avec ses contrôles.
- [x] Aucun moteur, champ de force ni paramètre de solvatation par défaut.

**Tests** : U (`test_structure`, `test_checks`, `test_config`), K (`test_interfaces`), architecture (`test_architecture`).

---

## Étape 1 — Modèle de données enrichi

**Objectif** : que `core` puisse représenter tout ce dont les étapes suivantes ont besoin, sans connaître de champ de force.

**Modules** : `core`, `residues`, `io.sequence`.

**Livrables**

- `Peptide` immuable (opérations qui renvoient un nouveau peptide : mutation, changement d'extrémité), chacune journalisée par l'appelant.
- Résidu avec état de protonation optionnel (`protonation: str | None`, clé interne neutre comme `his_delta`, `his_epsilon`, `his_protonated`).
- Modifications de chaînes latérales et liaisons non linéaires (pont disulfure, cyclisation tête-queue) déclarées, même si leur prise en charge complète vient plus tard.
- Notation des non standards dans une séquence (ex. `Ac-RG[Aib]S-NH2`), avec la règle de lecture documentée.

**Critères d'acceptation**

- [ ] Aucun nom de résidu de champ de force dans `core` (vérifié par un test qui cherche `NHE`, `NME`, `HIE`, `HID`, `HIP`, `CYX` dans `src/pepforge/core/`).
- [ ] Une mutation renvoie un nouveau `Peptide` ; l'original est inchangé.
- [ ] Les notations existantes (`Ac-RGDS-NH2`, `NH2-Arg-Gly-CO-NH2`...) donnent le même résultat qu'avant.

**Tests** : U (mutations, notation entre crochets, erreurs) ; V : aucune (pas encore de chimie quantitative).

**Revue** : relire la règle de lecture des séquences ; vérifier que les clés internes de protonation sont neutres.

---

## Étape 2 — Protonation et charge nette théorique

**Objectif** : choisir l'état de chaque groupe ionisable au pH de l'étude et calculer la charge nette théorique, référence de tous les contrôles de charge ultérieurs.

**Modules** : `residues/data/` (nouveau fichier de groupes ionisables), `residues`, `checks`, `pipeline` (étape « normalisation »).

**Prérequis** : choisir la **source** des pKa de référence (valeurs de résidus isolés publiées, ou outil de prédiction comme PROPKA pour une structure). À décider avec Firas ; aucune valeur n'est saisie sans sa source.

**Livrables**

- `residues/data/ionizable_groups.yaml` : groupes ionisables (amine N-terminale, carboxyle C-terminal, Asp, Glu, His, Cys, Tyr, Lys, Arg), charge de chaque état, pKa de référence et sa source.
- Fonction `theoretical_net_charge(peptide, ph) -> float` et choix des états (`assign_protonation(peptide, ph) -> Peptide`).
- Étape `ProtonationStage` et contrôles associés.

**Critères d'acceptation**

- [ ] Une extrémité acétylée ou amidée ne porte pas de charge ; une extrémité libre en porte une selon le pH.
- [ ] Le pH vient de `PrepConfig.ph` ; s'il manque, `ConfigurationError`.
- [ ] Un pH proche d'un pKa (fenêtre définie et justifiée) produit un `WARN` : l'état choisi doit être relu.
- [ ] Le choix de chaque état est journalisé avec la source du pKa.

**Tests** : U (chaque groupe, cas limites de pH) ; I (`SequenceInputStage` puis `ProtonationStage`) ; V : charges nettes des peptides de référence calculées à la main et consignées dans `tests/reference/charges.yaml` avec le raisonnement.

**Revue** : Firas valide la source des pKa et les charges attendues des peptides de référence.

---

## Étape 3 — Lecture et écriture PDB / mmCIF, inventaire

**Objectif** : importer une structure (expérimentale ou issue de docking) et en dresser l'inventaire, sans la modifier.

**Modules** : `io`, `core.structure`.

**Prérequis** : décider de lire les formats nous-mêmes ou via une bibliothèque (gemmi, Biopython...). Recommandation : une bibliothèque éprouvée pour mmCIF, encapsulée dans `io` pour que le reste du code ne voie que `Structure`.

**Livrables** : `io.read_structure(path) -> Structure`, `io.write_pdb(structure, path)`, `inventory(structure) -> StructureInventory` (chaînes, résidus, hétéroatomes, eaux, ions, conformations alternatives, hydrogènes, trous de numérotation).

**Critères d'acceptation**

- [ ] Lecture puis écriture puis relecture sans perte (coordonnées, noms, chaînes).
- [ ] L'inventaire signale tout ce qui demandera une décision (conformations alternatives, résidus manquants, ligands).
- [ ] Aucune correction silencieuse à la lecture.
- [ ] Correspondance `Structure` → `Peptide` quand la structure est un peptide reconnu.

**Tests** : U (petits fichiers PDB et mmCIF écrits à la main dans `tests/data/`) ; I (lecture d'une sortie de docking type) ; V : inventaire d'une structure publique d'un récepteur ciblé, comparé à la description de la PDB.

**Revue** : vérifier la gestion des conformations alternatives et des numéros d'insertion.

---

## Étape 4 — Contrôles de structure

**Objectif** : disposer de tous les contrôles de la section « Structure » de la checklist 04, applicables à n'importe quelle `Structure`.

**Modules** : `checks`, `residues/data/` (atomes attendus par résidu).

**Livrables** : contrôles « résidus attendus », « atomes lourds manquants », « hydrogènes présents », « conformations alternatives résolues », « extrémités conformes au peptide », « stéréochimie des carbones alpha », « chevauchements atomiques ». Données des atomes attendus tirées d'une source publiée (ex. dictionnaire des composants chimiques de la PDB), avec la version.

**Critères d'acceptation**

- [ ] Chaque item automatisable de la section « Structure » de la checklist 04 a un contrôle et l'indique dans `checklist_item`.
- [ ] Chaque contrôle a un test `PASS` et un test `FAIL`.
- [ ] Les seuils (distance de chevauchement...) viennent d'un fichier de données sourcé, pas du code.

**Tests** : U (structures construites à la main, correctes et volontairement défectueuses) ; V : structures publiques connues.

**Revue** : Firas relit les seuils et leurs sources.

---

## Étape 5 — Construction et correction 3D 🔒

**Décision requise** : outil de construction (ex. `tleap` d'AmberTools, constructeur Python, autre). Le choix n'affecte que `build`.

**Objectif** : produire une `Structure` à partir d'un `Peptide` (conformation étendue ou hélicoïdale) et corriger une structure importée (atomes manquants, hydrogènes selon la protonation de l'étape 2).

**Modules** : `build`, `tools` (si l'outil est externe).

**Livrables** : un `StructureBuilder`, des `StructureModifier` (sélection de chaîne, résolution des conformations alternatives, ajout d'atomes manquants, ajout des hydrogènes, mutation), étapes du pipeline correspondantes.

**Critères d'acceptation**

- [ ] `assert_builder_contract` passe.
- [ ] Tous les contrôles de l'étape 4 passent sur les peptides de référence construits.
- [ ] Extrémités conservées (vérifié par contrôle et par test).
- [ ] Chaque modification d'une structure importée est journalisée.

**Tests** : K, U, I (séquence → protonation → construction) ; V : géométrie (longueurs de liaison, angles dièdres de la conformation demandée) comparée à des valeurs de référence sourcées ; tests marqués `external` si un outil externe est utilisé.

**Revue** : inspection visuelle d'au moins une structure par peptide de référence (PyMOL ou VMD), consignée dans la PR.

---

## Étape 6 — Adaptateur CHARMM36m

**Champ de force retenu : CHARMM36m + CGenFF**, confirmé par Firas le 26 novembre 2026, cohérent avec NAMD. Modèle d'eau : TIP3P modifié CHARMM. Les fichiers à produire (`.rtf`, `.prm`, `.str`, `.psf`) et les deux pièges à éviter sont décrits en [ARCHITECTURE.md, section 11.3](ARCHITECTURE.md#113-les-fichiers-charmm-et-namd-à-produire).

**Objectif** : traduire un `Peptide` protoné dans la nomenclature du champ de force et produire un `ParameterSet` complet et sourcé.

**Modules** : `forcefields`, `forcefields/data/`.

**Livrables** : un adaptateur, ses tables de correspondance (résidus, extrémités, tautomères) avec leur source, la liste des modèles d'eau compatibles avec la publication qui l'établit.

**Critères d'acceptation**

- [ ] `assert_forcefield_contract` passe sur tous les peptides de référence.
- [ ] `-CO-NH2` donne `NHE` et `-CO-NH-CH3` donne `NME` en AMBER (ou l'équivalent documenté dans un autre champ de force) ; test dédié.
- [ ] Un résidu sans paramètres lève `MissingParameterError`.
- [ ] Un modèle d'eau incompatible lève `ConfigurationError`.

**Tests** : K, U (chaque correspondance) ; V : règle 7 de `CLAUDE.md`, charge nette et absence de paramètre manquant pour chaque résidu pris en charge.

**Revue** : Firas vérifie les correspondances de nomenclature et les citations.

---

## Étape 7 — Topologie pour NAMD

**Moteur retenu : NAMD**, d'après le projet de thèse (NAMD pour l'affinement des modèles, la dynamique moléculaire et la minimisation d'énergie, VMD pour l'outillage). Fichiers visés : `.psf` et `.pdb` via `psfgen`, plus les `.prm`, `.rtf` et `.str` du champ de force. Rappel : un `.str` CHARMM est un *stream file* de topologie et de paramètres, pas un fichier de coordonnées.

**Objectif** : écrire topologie, paramètres et coordonnées du soluté dans le vide.

**Modules** : `engines`, `tools`.

**Livrables** : un `EngineWriter` (`write_topology`), les motifs d'échec de l'outil utilisé (`failure_patterns`), les contrôles de l'étape 5 du pipeline (charge nette relue, nombre d'atomes, résidus terminaux).

**Critères d'acceptation**

- [ ] `assert_engine_writer_contract` passe.
- [ ] Charge nette **relue dans les fichiers produits** égale à la charge théorique de l'étape 2.
- [ ] Toute ligne d'avertissement de l'outil fait échouer l'étape, sauf exception documentée.
- [ ] Versions des outils dans `EngineFiles.tool_versions`.

**Tests** : K, U (génération des scripts d'entrée de l'outil), I (séquence → topologie, marqué `external`) ; V : énergie en point unique d'un peptide de référence comparée à un calcul fait à la main avec l'outil, consigné dans `tests/reference/`.

**Revue** : relire les fichiers d'entrée générés (scripts tleap, gmx ou psfgen).

---

## Étape 8 — Solvatation et neutralisation

**Objectif** : boîte, eau, contre-ions et sel selon `SolvationSpec`.

**Modules** : `engines` (`solvate`), `checks`, `config`.

**Critères d'acceptation**

- [ ] `SolvationSpec.missing_decisions()` vide, sinon `ConfigurationError` avant tout calcul.
- [ ] Système neutre ; concentration en sel conforme ; marge de boîte respectée ; aucun chevauchement (contrôles de la section « Système solvaté » de la checklist 04).
- [ ] Les choix (marge, sel) sont justifiés dans `PrepConfig.notes` et reportés dans le rapport.

**Tests** : U, I (marqués `external`) ; V : nombre d'ions calculé à la main pour un volume donné.

---

## Étape 9 — Protocole de simulation

**Objectif** : scripts de minimisation, d'équilibration (NVT puis NPT) et de production.

**Modules** : `config` (nouveau `ProtocolSpec`, sans valeur par défaut), `engines` (`write_protocol`), `engines/data/` (modèles de scripts commentés).

**Critères d'acceptation**

- [ ] Aucune valeur numérique dans le code : températures, durées, pas de temps, rayons de coupure viennent de `ProtocolSpec`, chacun avec sa justification.
- [ ] Contrôles de la section « Protocole de simulation » de la checklist 04 (pas d'intégration cohérent avec les contraintes sur les liaisons hydrogène, électrostatique longue portée déclarée, plusieurs répliques).
- [ ] Les scripts produits se terminent par leurs propres contrôles (règle 5 de `CLAUDE.md`).

**Tests** : U (rendu des modèles), V : exécution courte de chaque script sur un peptide de référence (marqué `external`, hors CI si trop long).

---

## Étape 10 — Rapport reproductible et peptide libre de bout en bout

**Objectif** : la commande `pepforge prepare config.yaml` enchaîne toutes les étapes pour un **peptide libre en solution** et produit un rapport rejouable.

**Modules** : `report`, `pipeline`, `cli`.

**Critères d'acceptation**

- [ ] Le rapport contient : configuration complète, journal, résultats de tous les contrôles, versions des outils, objectif de simulation, points à vérifier manuellement.
- [ ] Relancer avec le YAML du rapport donne des fichiers identiques (hors horodatages).
- [ ] Un échec à n'importe quelle étape produit un rapport partiel lisible.

**Tests** : I de bout en bout sur `Ac-RGDS-NH2` **et sur son analogue aza-β³ en une position** (dès que l'étape 11C est passée) ; test de reproductibilité (deux exécutions, fichiers comparés) ; V : courte simulation lancée et relue par Firas.

---

## Fil aza-β³ — étapes 11A à 11E (fonction centrale)

Ce fil réalise la substitution d'un acide aminé par un résidu aza-β³ et la paramétrisation qui va avec. Il démarre dès l'étape 4 passée, en parallèle des étapes 5 à 10.

> **Préalable, avant toute ligne de code : bibliographie.** Déterminer la variante exacte visée (aza-β³ `-NH-N(R)-CH2-CO-` ou azapeptide α `-NH-N(R)-CO-`), et chercher s'il existe des paramètres publiés pour des aza-peptides ou des β-peptides réutilisables, à quel niveau de théorie, et avec quelle validation conformationnelle. Livrable : une note dans `docs/` avec les références et la voie retenue, relue par Firas. Aucun paramètre n'est écrit dans le dépôt sans sa source ou sa procédure de dérivation complète.

### Étape 11A — Fragment chimique et template de résidu

**Modules** : `core` (`ResidueTemplate`), `io` (SMILES, MOL/SDF), `checks`.

**Décision requise** : bibliothèque chimique (RDKit ou autre).

**Livrables** : type `ResidueTemplate` (fragment, atomes d'ancrage côté N et côté C, atomes retirés au greffage, charge attendue), lecture d'un fragment depuis SMILES ou SDF, contrôles de valence, de stéréochimie et de protonation du fragment.

**Critères d'acceptation**

- [ ] Les points d'ancrage sont déclarés explicitement, jamais devinés.
- [ ] Un fragment de valence incohérente lève `ChemistryError`.
- [ ] Le template ne contient aucun paramètre de champ de force : c'est de la chimie, pas de la paramétrisation.

**Tests** : U (fragments corrects et défectueux) ; K.

### Étape 11B — Greffage dans le peptide

**Modules** : `build` (`MutateResidue`), `checks`, `pipeline`.

**Livrables** : substitution du résidu en position i par un `ResidueTemplate`, reconstruction des liaisons de jonction, minimisation locale du voisinage, journalisation complète.

**Critères d'acceptation**

- [ ] Le squelette reste continu après greffage (contrôle de connectivité, et non un comptage d'atomes supposé).
- [ ] Le reste de la structure est inchangé hors de la zone minimisée, et les extrémités du peptide sont conservées.
- [ ] La longueur et les angles de la liaison N–N tombent dans une plage issue d'une source citée.
- [ ] Les contrôles de stéréochimie attendent un Cα chiral sur les résidus α et aucun centre stéréogène sur le résidu aza.
- [ ] Chaque substitution est journalisée : position, résidu retiré, template ajouté, paramètres de minimisation.

**Tests** : U, K, I (séquence → protonation → construction → substitution) ; V : inspection visuelle de chaque peptide de référence substitué, consignée dans la PR.

### Étape 11C — Dérivation des paramètres d'un résidu aza-β³ 🔒

**Décisions requises** : voie de paramétrisation (AMBER + GAFF2 ou CHARMM + CGenFF), méthode de charges et niveau de théorie, logiciel de calcul quantique disponible.

**Modules** : `forcefields` (nouveau sous-module `derivation`), `tools`, `report`.

**Livrables** : chaîne fragment coiffé (type `Ace-X-Nme`) → calcul quantique → charges → `parmchk2` ou CGenFF → fichier de paramètres réutilisable, avec niveau de confiance et procédure complète.

**Critères d'acceptation**

- [ ] La charge nette du résidu dérivé est entière, et les atomes de jonction restent compatibles avec les résidus voisins du champ de force protéique.
- [ ] Tout paramètre obtenu par analogie est signalé avec son niveau de confiance (pénalités CGenFF, avertissements `parmchk2`).
- [ ] Aucune torsion manquante autour de la liaison N–N.
- [ ] La dérivation est rejouable à l'identique depuis le seul fichier de configuration (règle 8), et le rapport consigne méthode de charges, niveau de théorie, pénalités et versions des outils.

**Tests** : U, K ; V : règle 7 (charge nette, aucun paramètre manquant) pour chaque résidu dérivé.

### Étape 11D — Validation conformationnelle des paramètres dérivés 🔒

**Objectif** : sans cette étape, les conformations produites par la simulation ne sont pas interprétables, et les paramètres ne sont pas défendables dans une thèse.

**Critères d'acceptation**

- [ ] Balayages dièdres relaxés en calcul quantique pour φ, ψ et la torsion N–N d'un résidu aza-β³ de référence.
- [ ] Courbes champ de force et courbes quantiques comparées, écart maximal consigné, seuil d'acceptation fixé et justifié avec Firas.
- [ ] Courbes et données conservées dans `tests/reference/`.
- [ ] Si l'écart dépasse le seuil : ajustement des torsions documenté, ou point signalé comme limite connue dans le rapport.

**Tests** : V, avec relecture scientifique obligatoire de Firas.

### Étape 11E — Bibliothèque complète et scan de positions

**Objectif** : « fonctionner sur tous les acides aminés », c'est-à-dire une entrée par chaîne latérale **et par état de protonation**, soit une trentaine d'entrées, puis la substitution en lot.

**Modules** : `residues/data/`, `pipeline` (mode lot), `cli`.

**Critères d'acceptation**

- [ ] Une entrée de bibliothèque par chaîne latérale et par état de protonation pertinent, chacune avec sa procédure de dérivation et sa version.
- [ ] Les cas particuliers sont traités explicitement, et non ignorés : glycine (pas de chaîne latérale), proline (azote cyclique), cystéine (pont disulfure), chaînes latérales ionisables.
- [ ] Le mode lot prépare un système par position pour i = 1..n, chacun avec son rapport.
- [ ] Chaque entrée ajoutée passe le test de charge nette et d'absence de paramètre manquant (règle 7) avant d'entrer dans la bibliothèque.

**Tests** : U par entrée ; I sur un scan complet d'un peptide court ; V : relecture par Firas des entrées, chaîne latérale par chaîne latérale.

---

## Étape 11F — Chélateurs et radiomarquage 🔒

**Périmètre retenu : les deux cas**, décidé par Firas le 26 novembre 2026 — le peptide portant son chélateur **sans métal**, et le peptide **marqué**. Le volet sans métal suit la mécanique de l'étape 11A. Le volet métal est d'une autre nature : voir [ARCHITECTURE.md, section 11.4](ARCHITECTURE.md#114-le-centre-métallique--un-problème-de-nature-différente).

**Décision encore requise** : l'approche de paramétrisation du centre métallique (non lié, lié, atomes factices cationiques) et le format retenu pour le système marqué. Elle se prend après la recherche bibliographique, pas avant.

**Contexte** : les peptides de la thèse portent un chélateur greffé à l'extrémité N-terminale ou C-terminale amide, voire sur une chaîne latérale : DOTA (DOTA(tBu)₃), N4 ((Boc)₄-N4), HYNIC (Boc-HYNIC). Radiométaux visés : gallium, rhénium, lutécium, technétium-99m.

**Modules** : `residues/data/` (chélateurs comme `ResidueTemplate`), `forcefields` (dérivation), `checks`.

**Critères d'acceptation, volet chélateur sans métal**

- [ ] Un chélateur est décrit comme un `ResidueTemplate` terminal ou latéral, avec ses points d'ancrage : même mécanique que l'étape 11A, aucun code spécifique.
- [ ] Greffage en N-terminal, en C-terminal amide et sur une chaîne latérale, chacun testé.
- [ ] Forme protégée (tBu, Boc) et forme déprotégée distinguées explicitement : seule la forme présente dans le peptide final est simulée.
- [ ] Charge nette du chélateur conforme à son état de protonation au pH de l'étude.

**Critères d'acceptation, volet radiométal (seulement si la décision de périmètre l'inclut)**

- [ ] Approche de paramétrisation du centre métallique choisie et justifiée avec sa source parmi les trois familles de la section 11.4, avec ses limites écrites dans le rapport.
- [ ] Aucun mélange entre champ de force protéique CHARMM et paramètres métalliques AMBER dans un même système : un contrôle le refuse (règle 4 du projet).
- [ ] Pour un cœur oxo de technétium ou de rhénium, l'hypothèse structurale retenue est énoncée explicitement dans le rapport, puisque la nature du complexe formé au radiomarquage n'est pas toujours établie.
- [ ] Géométrie de coordination contrôlée après construction.
- [ ] Charge totale du complexe vérifiée contre la charge théorique attendue.
- [ ] Limites de l'approche signalées dans le rapport de chaque préparation concernée.

**Tests** : U, K ; V : comparaison avec une structure de complexe publiée pour le même chélateur et le même métal, relue par Firas.

---

## Étape 12 — Complexe peptide–récepteur

**Objectif** : préparer le récepteur (DCC, TrkC, UNC5, récepteur de la nétrine-1), placer le peptide issu du docking, paramétrer l'ensemble. **Priorité relevée** : le projet de thèse enchaîne docking puis dynamique moléculaire pour évaluer affinité et spécificité, donc le complexe est un cas principal et non une extension.

**Modules** : `core` (plusieurs chaînes), `io`, `build`, `pipeline` (variante de chaîne selon `PrepConfig.goal`).

**Critères d'acceptation**

- [ ] `goal = peptide_receptor_complex` sans récepteur : `ConfigurationError`.
- [ ] Le récepteur passe les mêmes contrôles de structure que le peptide.
- [ ] La taille de boîte et le protocole tiennent compte de l'objectif ; les deux scénarios ne partagent pas de valeurs implicites.

**Tests** : I sur un complexe de référence ; V : relecture de l'interface et des contacts par Firas.

---

## Étape 13 — Second moteur, second champ de force, conversion

**Objectif** : prouver que l'architecture tient : ajouter un moteur ou un champ de force sans modifier les modules existants, et convertir entre moteurs.

**Critères d'acceptation**

- [ ] L'ajout ne modifie ni `core`, ni `build`, ni les autres adaptateurs (vérifiable dans le diff de la PR).
- [ ] Conversion traitée comme une étape distincte : nombre d'atomes, charge nette et énergie en point unique comparés avant et après.

**Tests** : K pour la nouvelle implémentation ; V : comparaison croisée des énergies entre moteurs.

---

## Étape 14 — Analyses

**Objectif** : analyses de la section 7 de la spécification (RMSD, RMSF, rayon de giration, structure secondaire, liaisons hydrogène, SASA, contacts d'interface, regroupement).

**Architecture** : nouveau module `analysis`, hors de la chaîne de préparation. Il consomme des trajectoires et le rapport de préparation, et ne dépend que de `core`, `config`, `io` et `report`. À déclarer dans la matrice de dépendances.

**Tests** : U sur des trajectoires synthétiques ; V : comparaison avec un outil d'analyse publié sur une même trajectoire.
