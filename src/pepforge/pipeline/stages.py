"""Étapes concrètes du pipeline, ajoutées au fil du plan de construction."""

from __future__ import annotations

from dataclasses import dataclass

from pepforge.core.peptide import Peptide
from pepforge.io.sequence import parse_sequence
from pepforge.pipeline.stage import StageContext
from pepforge.residues.library import ResidueLibrary


@dataclass
class SequenceInputStage:
    """Étape 1 (mode séquence) : texte → `Peptide`."""

    library: ResidueLibrary | None = None
    name: str = "input.sequence"

    def run(self, data: str, context: StageContext) -> Peptide:
        peptide = parse_sequence(data, self.library)
        context.journal.record(
            "parse_sequence",
            "lecture de la séquence d'entrée",
            input=data,
            n_terminus=peptide.n_terminus.key,
            c_terminus=peptide.c_terminus.key,
            length=len(peptide),
        )
        return peptide
