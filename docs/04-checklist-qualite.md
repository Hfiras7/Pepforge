# PepForge — Contrôles qualité avant simulation

Liste à parcourir pour chaque système préparé. Un point non vérifié est un point en échec.

## Structure

- [ ] Nombre de résidus conforme à la séquence voulue.
- [ ] Extrémités N- et C-terminales conformes à la notation d'entrée. Pour une extrémité `-CO-NH2` en AMBER, le résidu terminal doit être `NHE` et non `NME`.
- [ ] Aucun atome ni résidu manquant après reconstruction.
- [ ] Liaisons hors squelette (ponts disulfure, cyclisation) et modifications de chaîne latérale : chimiquement possibles, vérifiables, et présentes dans la topologie produite.
- [ ] Hydrogènes présents et cohérents avec le pH retenu.
- [ ] Stéréochimie des centres chiraux préservée.
- [ ] Pour une structure importée : conformations alternatives résolues, solvants et ligands traités de façon explicite.

## Paramétrisation

- [ ] Journal de construction relu : aucun message de paramètre manquant.
- [ ] Charge nette calculée identique à la charge théorique attendue au pH choisi.
- [ ] Types atomiques attribués pour tous les résidus, non standards inclus.
- [ ] Pour les fragments personnalisés : procédure de dérivation des charges documentée, paramètres obtenus par analogie signalés.
- [ ] Champ de force et modèle d'eau compatibles entre eux et avec le moteur visé.

## Système solvaté

- [ ] Marge entre le soluté et le bord de la boîte suffisante pour éviter les interactions avec les images périodiques.
- [ ] Système électriquement neutre après ajout des contre-ions.
- [ ] Concentration en sel conforme à celle visée.
- [ ] Aucun chevauchement atomique après solvatation.

## Protocole de simulation

- [ ] Minimisation convergée, énergie finale sans valeur aberrante.
- [ ] Équilibration en deux temps au moins : volume constant puis pression constante.
- [ ] Température, pression, densité et énergie stables en fin d'équilibration.
- [ ] Pas d'intégration cohérent avec le traitement des liaisons impliquant les hydrogènes.
- [ ] Traitement des interactions électrostatiques à longue portée explicite, rayons de coupure déclarés.
- [ ] Plusieurs répliques indépendantes prévues, avec des vitesses initiales différentes.

## Conversion entre moteurs

- [ ] Conversion traitée comme une étape à part, jamais implicite.
- [ ] Après conversion : nombre d'atomes, charge nette et énergie en point unique comparés entre les deux formats.

## Traçabilité

- [ ] Versions de tous les outils externes enregistrées.
- [ ] Fichier de configuration permettant de relancer la préparation à l'identique.
- [ ] Journal complet des modifications apportées à la structure de départ.
- [ ] Objectif de la simulation consigné : peptide libre en solution, ou complexe peptide–récepteur.
