"""Défauts injectés un à un dans le module 3 : chacun doit faire échouer au moins un test.

    python validation/module03_mutations.py

Chaque défaut désactive une garde (refus, comparaison, plafond) ; le fichier est
restauré après chaque essai. Lancé avec PYTHONDONTWRITEBYTECODE=1 pour éviter
un bytecode périmé (voir le bilan du module 1).
"""

import os
import pathlib
import subprocess
import sys

M = [
    ("src/pepforge/io/structure_formats.py", "if not _FLOAT.fullmatch(value):", "if False:"),
    (
        "src/pepforge/io/structure_formats.py",
        '    if value in ("", "?", "."):\n        return None',
        '    if value in ("", "?", "."):\n        return 1.0',
    ),
    (
        "src/pepforge/io/structure_formats.py",
        "        if item.key in expected:\n            raise",
        "        if False:\n            raise",
    ),
    ("src/pepforge/io/structure_formats.py", "    if lost:\n", "    if False:\n"),
    (
        "src/pepforge/io/structure_formats.py",
        "                    if source is None or key in seen:",
        "                    if source is None:",
    ),
    ("src/pepforge/io/structure_formats.py", "abs_tol=_COORDINATE_TOLERANCE", "abs_tol=1.0"),
    (
        "src/pepforge/io/structure_formats.py",
        '        if not started:\n            started = "RES C SSSEQI" in body\n            continue',
        "        if not started:\n            started = True\n            continue",
    ),
    (
        "src/pepforge/io/structure_formats.py",
        'if get("polymer_flag") != "Y" or get("occupancy_flag") != "1":',
        'if get("polymer_flag") != "Y":',
    ),
    (
        "src/pepforge/io/structure_formats.py",
        '    if in_model:\n        raise StructureFormatError(f"{origin} : dernier',
        '    if False:\n        raise StructureFormatError(f"{origin} : dernier',
    ),
    (
        "src/pepforge/io/structure.py",
        "        if len(numbers) > 1:\n            raise",
        "        if len(numbers) > 99:\n            raise",
    ),
    (
        "src/pepforge/io/structure.py",
        "if parsed.declared_model_count is not None and parsed.declared_model_count != len(numbers):",
        "if False:",
    ),
    (
        "src/pepforge/io/structure.py",
        "        if value is not None and value != by_content:",
        "        if False:",
    ),
    (
        "src/pepforge/io/structure.py",
        "    if structure_format is StructureFormat.PDB:\n        problems",
        "    if False:\n        problems",
    ),
    (
        "src/pepforge/io/structure.py",
        "        _compare_written(structure.atoms, written, structure_format, str(target))\n",
        "\n",
    ),
    (
        "src/pepforge/io/structure.py",
        "        os.replace(temporary, target)\n    except BaseException:\n        Path(temporary).unlink(missing_ok=True)",
        "        os.replace(temporary, target)\n    except BaseException:\n        pass",
    ),
    ("src/pepforge/io/retrieval.py", "        if total > max_bytes:", "        if False:"),
    (
        "src/pepforge/io/retrieval.py",
        "            if announced is not None and announced > limits.max_bytes:",
        "            if False:",
    ),
    ("src/pepforge/io/retrieval.py", '"text/html",', '"text/nothing",'),
    (
        "src/pepforge/io/retrieval.py",
        'allowed = {"https", "http"} if limits.allow_http else {"https"}\n    scheme',
        'allowed = {"https", "http"}\n    scheme',
    ),
    ("src/pepforge/io/retrieval.py", "    if parts.username or parts.password:", "    if False:"),
    (
        "src/pepforge/io/retrieval.py",
        "    if announced is not None and len(raw) != announced:",
        "    if False:",
    ),
    (
        "src/pepforge/io/retrieval.py",
        "    if len(content) > max_bytes or decompressor.unconsumed_tail:",
        "    if False:",
    ),
    (
        "src/pepforge/io/retrieval.py",
        'if scheme not in self.allowed or (original == "https" and scheme != "https"):',
        "if scheme not in self.allowed:",
    ),
    ("src/pepforge/io/retrieval.py", "            if status != 200:", "            if False:"),
    (
        "src/pepforge/io/retrieval.py",
        "sha256=hashlib.sha256(raw).hexdigest()",
        "sha256=hashlib.sha256(content).hexdigest()",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "        and name not in components.c_terminal_caps\n",
        "\n",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "    if missing:\n        raise StructureError(",
        "    if False:\n        raise StructureError(",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "        if after.number - before.number > 1 or after.number < before.number:\n            raise",
        "        if False:\n            raise",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "            if first is None or second is None:\n                raise UnsupportedFeatureError",
        "            if False:\n                raise UnsupportedFeatureError",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "    if c_key is None and last_names & components.free_c_terminus_atoms:",
        "    if c_key is None and False:",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "    if inv.model_count and inv.model_count > 1:",
        "    if False:",
    ),
    (
        "src/pepforge/io/structure_inventory.py",
        "        if site.residue.chain_id == chain_id and site.residue_names:",
        "        if False:",
    ),
    (
        "src/pepforge/core/structure.py",
        "            if not math.isfinite(value):",
        "            if False:",
    ),
]
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
caught = 0
for i, (path, old, new) in enumerate(M, 1):
    p = pathlib.Path(path)
    src = p.read_text()
    assert src.count(old) == 1, (i, path, old)
    p.write_text(src.replace(old, new))
    try:
        r = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "-x",
                "-p",
                "no:cacheprovider",
                "tests/test_structure_io.py",
                "tests/test_retrieval.py",
                "tests/test_structure_inventory.py",
                "tests/test_cli.py",
            ],
            capture_output=True,
            text=True,
            env=env,
        )
        ok = r.returncode != 0
        failing = [line for line in r.stdout.splitlines() if line.startswith("FAILED")][:1]
        caught += ok
        print(
            f"{i:2d} {'DÉTECTÉ ' if ok else 'NON DÉTECTÉ'} {path.split('/')[-1]}: {new.strip()[:50]!r} -> {failing[0][7:90] if failing else ''}"
        )
    finally:
        p.write_text(src)
print(f"{caught}/{len(M)} défauts détectés")
