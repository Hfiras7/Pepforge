# PepForge — Spécification fonctionnelle

## 1. Identité de l'outil

- **Nom** : PepForge
- **Sous-titre** : *Universal Peptide Design, Engineering and Simulation Platform*
- **Positionnement** : plateforme de conception, de modification, de construction 3D et de préparation des peptides pour la dynamique moléculaire.
- **À vérifier avant publication** : disponibilité du nom comme nom de logiciel et comme nom de domaine. Alternatives écartées mais conservées en réserve : PepEngine, PepCraft, PepParam, PepStruct.

## 2. Modes d'entrée

PepForge accepte trois points de départ. Chaque mode doit déboucher sur la même représentation interne du peptide.

### 2.1 Séquence peptidique

- Saisie en code à une lettre ou à trois lettres.
- Reconnaissance d'une notation complète incluant les extrémités, par exemple `NH2-...-CO-NH2`.
- Déclaration explicite des extrémités N- et C-terminales (libre, acétylée, amidée, etc.).
- Insertion d'acides aminés non standards et de modifications de chaînes latérales depuis une bibliothèque de résidus.
- Calcul automatique de la composition, de la longueur et de la charge nette théorique au pH choisi.
- Construction 3D initiale (conformation étendue ou hélicoïdale au choix) avec génération des hydrogènes.

### 2.2 Import d'une structure 3D

- Lecture de fichiers PDB, mmCIF et, au besoin, de sorties de docking moléculaire.
- Inventaire de la structure : chaînes, résidus manquants, atomes manquants, hydrogènes, hétéroatomes, molécules d'eau et ions.
- Nettoyage contrôlé : sélection de chaînes, suppression ou conservation des solvants et ligands, choix des conformations alternatives.
- Reconstruction des atomes manquants et ajout des hydrogènes selon le pH.
- Modification *a posteriori* : mutation de résidus, changement d'extrémités, introduction de modifications.
- Conservation d'un journal de toutes les modifications appliquées à la structure importée.

### 2.3 Structure chimique personnalisée

Mode destiné aux résidus et aux modifications absents des champs de force standards.

- Entrée par SMILES, MOL/SDF ou dessin de la structure.
- Vérification de la valence, de la stéréochimie et de l'état de protonation.
- Définition des points d'ancrage permettant l'intégration du fragment dans la chaîne peptidique.
- Attribution des types atomiques et dérivation des charges par la voie appropriée au champ de force cible (GAFF2 via `antechamber`/`parmchk2` pour AMBER, CGenFF pour CHARMM).
- Signalement explicite de tout paramètre manquant ou obtenu par analogie, avec la mention du niveau de confiance.
- Production d'un fichier de paramètres de résidu réutilisable (bibliothèque locale versionnée).

## 3. Sorties attendues

| Famille | Fichiers | Moteur |
|---|---|---|
| Structure | `.pdb`, `.mmcif` | tous |
| Topologie et paramètres AMBER | `.prmtop`, `.inpcrd` | AMBER, GROMACS après conversion |
| Topologie GROMACS | `.top`, `.itp`, `.gro` | GROMACS |
| Topologie et paramètres CHARMM | `.psf`, `.prm`/`.par`, `.str`, `.rtf` | NAMD, CHARMM |
| Scripts de préparation | solvatation, neutralisation, minimisation, équilibration, production | selon le moteur |
| Rapport | journal des étapes, charge nette, paramètres manquants, versions des outils | — |

> Note sur la terminologie employée dans la discussion initiale : un fichier `.str` CHARMM est un *stream file* (topologie et paramètres d'un résidu), et non un fichier de coordonnées. Les coordonnées se trouvent dans le `.pdb` ou le `.crd`. Cette distinction doit être respectée dans l'interface et la documentation.

## 4. Chaîne de traitement

1. Entrée (séquence, structure 3D ou structure personnalisée).
2. Normalisation et validation chimique.
3. Choix du couple champ de force / moteur.
4. Attribution des paramètres, résidus personnalisés inclus.
5. Construction ou correction de la structure 3D.
6. Génération de la topologie et des paramètres.
7. Contrôles automatiques (voir `04-checklist-qualite.md`).
8. Solvatation, ajout des ions, neutralisation.
9. Génération des scripts de minimisation, d'équilibration et de production.
10. Rapport de préparation reproductible.

## 5. Champs de force et moteurs à prendre en charge

- **AMBER** : ff14SB, ff19SB, GAFF2 pour les fragments non standards. Modèles d'eau TIP3P, OPC.
- **CHARMM** : CHARMM36m, CGenFF. Modèle d'eau TIP3P modifié CHARMM.
- **GROMOS** : prise en charge optionnelle, à prioriser plus tard.
- **Moteurs** : AMBER (`pmemd`), GROMACS, NAMD. Les conversions entre formats sont des étapes explicites et vérifiées, jamais implicites.

## 6. Deux objectifs de simulation distincts

1. **Peptide libre en solution** : étude de la stabilité conformationnelle, des structures secondaires et du repliement.
2. **Complexe peptide–récepteur** : préparation du récepteur (DCC, TrkC, UNC5), placement du peptide issu du docking, paramétrisation de l'ensemble, étude de la stabilité du complexe et des contacts d'interface.

L'interface doit demander lequel des deux est visé, car la préparation, la taille de la boîte, la durée de simulation et les analyses diffèrent.

## 7. Analyses prévues

RMSD, RMSF, rayon de giration, structure secondaire au cours du temps, liaisons hydrogène, surface accessible au solvant, contacts d'interface et énergie d'interaction pour les complexes, regroupement conformationnel.
