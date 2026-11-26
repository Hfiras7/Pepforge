# Prompt — Génération des fichiers d'entrée de dynamique moléculaire

> Prompt réutilisable. Remplacer les champs entre crochets puis l'envoyer à l'assistant.

---

## Prompt

Tu agis comme expert en préparation de systèmes de dynamique moléculaire pour PepForge.

**Entrée disponible**
- Type d'entrée : [séquence seule / fichier PDB / PDB issu du docking / structure chimique personnalisée]
- Peptide : [séquence avec extrémités, par exemple `NH2-...-CO-NH2`]
- Modifications non standards : [aucune / liste]
- pH de référence : [7,4 par défaut]

**Cible logicielle**
- Moteur : [AMBER / GROMACS / NAMD]
- Champ de force : [ff14SB / ff19SB / CHARMM36m / à déterminer selon le moteur]
- Modèle d'eau et ions : [TIP3P + NaCl 0,15 M par défaut]

**Objectif**
- [peptide libre en solution / complexe peptide–récepteur / les deux]
- Récepteur, le cas échéant : [DCC / TrkC / UNC5 + code PDB ou fichier]
- Durée de production visée : [par exemple 500 ns, 3 répliques]

**Ce que je veux en retour, dans cet ordre**

1. **Analyse de la séquence** : longueur, composition en résidus ionisables, états de protonation retenus au pH indiqué, charge nette théorique attendue.
2. **Points de vigilance** : traitement des extrémités, résidus nécessitant un paramétrage particulier, pièges de nomenclature propres au champ de force choisi.
3. **Script de construction complet et commenté** (`tleap` pour AMBER, `pdb2gmx` ou équivalent pour GROMACS, `psfgen` pour NAMD/CHARMM), incluant les commandes de vérification de structure et de charge avant la sauvegarde des fichiers.
4. **Script de solvatation et de neutralisation** : type de boîte, marge au solvant, nombre d'ions, concentration en sel.
5. **Protocole de minimisation, d'équilibration et de production** : fichiers d'entrée du moteur choisi, avec pour chaque étape l'ensemble thermodynamique, les contraintes, le pas d'intégration, le thermostat, le barostat, le traitement des interactions à longue portée et la durée.
6. **Contrôles à effectuer** : ce qu'il faut lire dans le journal de construction, et les signes qu'une étape a échoué.
7. **Analyses de post-traitement** : commandes pour RMSD, RMSF, rayon de giration, structure secondaire, liaisons hydrogène, et pour un complexe, contacts d'interface.
8. **Versions des outils** supposées, et ce qui change avec une version différente.

**Contraintes**
- Rédige en français ; garde les noms de commandes, de fichiers et de résidus dans leur forme originale.
- N'invente aucun paramètre de champ de force. Pour tout résidu non standard, décris la procédure de dérivation et indique le niveau de confiance.
- Signale explicitement toute conversion de format entre moteurs comme une étape à part, avec sa vérification.
- Si une information manquante change le résultat, pose la question avant de produire le workflow.

---

## Exemple d'ossature attendue pour la construction AMBER

```tleap
source leaprc.protein.ff14SB

# Construction du peptide avec les extrémités voulues
pep = sequence { <residus avec residu terminal adapte> }

# Verifications avant sauvegarde
check pep
charge pep

# Fichiers de sortie
saveamberparm pep peptide.prmtop peptide.inpcrd
savepdb pep peptide.pdb

quit
```

Le journal de LEaP doit être relu pour confirmer que les deux extrémités ont été interprétées comme voulu, qu'aucun paramètre n'est manquant, et que la charge obtenue correspond à l'état de protonation choisi.
