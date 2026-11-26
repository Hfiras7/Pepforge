## Étape du plan

<!-- Numéro et titre de l'étape de docs/PLAN-CONSTRUCTION.md. -->

## Ce qui change

<!-- Résumé en quelques lignes. -->

## Critères d'acceptation de l'étape

<!-- Recopier les critères de l'étape et les cocher un à un. -->

- [ ] ...

## Modules touchés et dépendances

- Modules modifiés :
- [ ] La matrice de dépendances (`tests/test_architecture.py`) n'a pas changé, **ou** son changement est justifié ici et reporté dans `docs/ARCHITECTURE.md`.
- [ ] Les contrats et types d'échange décrits dans `docs/ARCHITECTURE.md` sont à jour.

## Tests ajoutés

- [ ] Unitaires :
- [ ] Contrat (`tests/contracts.py`) :
- [ ] Intégration :
- [ ] Validation scientifique (valeurs de référence et leur source) :

## Données et choix scientifiques

<!-- Toute nouvelle donnée (paramètre, pKa, seuil, correspondance de nomenclature) avec sa source. -->

- [ ] Aucun paramètre inventé ; chaque donnée a sa source ou est signalée « à déterminer ».
- [ ] Toute opération modifiant une structure est journalisée.

## Vérifications

- [ ] `ruff check . && ruff format --check .`
- [ ] `mypy`
- [ ] `pytest`
- [ ] Relecture scientifique par Firas (si l'étape le demande)
