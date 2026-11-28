# PepForge

*Universal Peptide Design, Engineering and Simulation Platform*

PepForge est une plateforme Python de conception, de modification, de construction 3D et de préparation de peptides pour la dynamique moléculaire. Elle part d'une séquence, d'une structure 3D (PDB, mmCIF, sortie de docking) ou d'une structure chimique personnalisée, et produit les fichiers d'entrée de simulation (topologie, paramètres, coordonnées) avec les scripts de préparation associés.

Le projet est développé dans le cadre d'une thèse de doctorat sur des peptides théranostiques ciblant les récepteurs à dépendance (DCC, TrkC, UNC5).

> **État actuel : squelette.** L'analyse des séquences et de leurs extrémités, la configuration, le journal et l'encapsulation des outils externes sont en place. La construction 3D, les adaptateurs de champ de force et la génération des fichiers moteur restent à écrire. Le moteur de simulation, le champ de force et le système simulé ne sont pas encore choisis et restent configurables.

## Installation

Prérequis : Python 3.11 ou plus, et Git.

```bash
git clone https://github.com/Hfiras7/Pepforge.git
cd Pepforge
python -m venv .venv
source .venv/bin/activate        # Windows : .venv\Scripts\activate
pip install -e ".[dev]"
```

## Utilisation

```bash
# Le mode de lecture se précise toujours : 1 = code à une lettre, 3 = code à trois lettres
pepforge parse --notation 1 "Ac-RGDS-NH2"
pepforge parse --notation 3 "NH2-Arg-Gly-Asp-Ser-CO-NH2"
pepforge convert --notation 1 "Ac-RGDS-NH2"      # → Ac-Arg-Gly-Asp-Ser-NH2
pepforge convert --notation 3 "Ac-Ala-NHMe"      # → Ac-A-NHMe
pepforge config-template > prep.yaml

# Importer une structure, sans la modifier, et afficher son inventaire
pepforge inspect structure.cif
pepforge inspect ensemble.pdb --model 1 --chain A     # modèle à choisir si plusieurs
pepforge inspect --url https://exemple.org/structure.cif
pepforge inspect structure.cif --chain A              # rapport d'extraction : jonctions C–N, composants, refus
pepforge inspect structure.cif --chain A --breaks split   # décision explicite, journalisée
```

## Développement

```bash
pytest                                  # tests
ruff check . && ruff format --check .   # style
mypy                                    # typage strict
```

### Avec VS Code et Claude Code

1. Ouvrir le dossier du dépôt dans VS Code.
2. Accepter l'installation des extensions recommandées (Python, Pylance, Ruff, Mypy, Claude Code).
3. Sélectionner l'interpréteur `.venv` (palette de commandes : *Python: Select Interpreter*).
4. Lancer Claude Code et se connecter avec son compte Claude. Il lit automatiquement [`CLAUDE.md`](CLAUDE.md), qui contient les règles du projet.

## Organisation

```
src/pepforge/
├── core/         représentation interne (Peptide, Structure)
├── config.py     configuration d'une préparation (sans valeur par défaut)
├── tools/        encapsulation des exécutables externes
├── checks/       contrôles qualité après chaque étape
├── build/        construction et correction 3D (interfaces)
├── io/           séquences, puis PDB, mmCIF, SMILES, MOL/SDF
├── residues/     bibliothèque de résidus et extrémités (données YAML sourcées)
├── forcefields/  adaptateurs de champ de force (interface commune)
├── engines/      fichiers d'entrée AMBER, GROMACS, NAMD (interface commune)
├── pipeline/     enchaînement des étapes et contrôles
├── report/       journal et rapport reproductible
└── cli/          commande `pepforge`
docs/             instructions, spécifications, ARCHITECTURE.md, PLAN-CONSTRUCTION.md
tests/            tests pytest
```

## Documentation

- [Instructions du projet](docs/00-instructions-projet.md)
- [Spécification fonctionnelle](docs/01-specification-fonctionnelle.md)
- [Prompt de génération des fichiers de DM](docs/02-prompt-generation-fichiers-md.md)
- [Prompt de développement logiciel](docs/03-prompt-developpement-logiciel.md)
- [Checklist qualité avant simulation](docs/04-checklist-qualite.md)
