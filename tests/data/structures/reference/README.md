# Structures de référence pour les tests d'extraction (module 3)

Trois fichiers réels, compressés (`gzip -n -9`) sans aucune modification de leur
contenu, pour que les tests d'extraction tournent sans réseau. Ils viennent de la
suite de tests de Biopython 1.88, elle-même copie d'entrées de la Protein Data Bank
(données PDB en CC0 1.0) :
`https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/PDB/<fichier>`.
Ce sont les mêmes octets que ceux de la validation (`validation/module03_structures.py`,
`validation/results/module03_structures.yaml`).

| Fichier | SHA-256 du fichier décompressé | Cas couvert |
|---|---|---|
| `4ZHL.cif.gz` | `b40655504ff4982cf8df49911019ada322d410a93d9bf66c65a9000e2531f60f` | numérotation chymotrypsine (saut ARG217 → GLY219 avec liaison C–N), codes d'insertion |
| `6WQA.cif.gz` | `484ff1689177aee680604631b9139199db8064f7de783e267f51cdfdf1769e6f` | protéine de fusion (sauts +793 et −887 avec liaison), rupture non déclarée ALA1043 → ASP1060, ligand non lié `ZMA` |
| `7DDO.pdb.gz` | `d6f4f7bacd3a8c8d21c4ec9d2543d74887a38fe469fd9840f8b9275f4039dd50` | glycanes `NAG` liés aux Asn (`LINK`), ion Zn coordonné, résidus non observés aux extrémités (chaîne C) |

Les tests vérifient l'empreinte de chaque fichier décompressé avant de s'en servir.
