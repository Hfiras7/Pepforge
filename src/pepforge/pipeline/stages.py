"""Étapes concrètes du pipeline, ajoutées au fil du plan de construction."""

from __future__ import annotations

from dataclasses import dataclass

from pepforge.core.peptide import Peptide
from pepforge.core.protonation import C_TERMINUS_SITE, N_TERMINUS_SITE, ProtonationResult
from pepforge.exceptions import ConfigurationError
from pepforge.io.sequence import SEQUENCE_NOTATION, SequenceNotation, parse_sequence
from pepforge.pipeline.stage import StageContext
from pepforge.residues.library import ResidueLibrary
from pepforge.residues.protonation import (
    IonizableGroupTable,
    assign_protonation,
    load_ionizable_groups,
)


@dataclass
class SequenceInputStage:
    """Étape 1 (mode séquence) : texte → `Peptide`, lu dans le code déclaré au départ.

    Le mode de lecture n'a pas de valeur par défaut : code à une lettre ou code
    à trois lettres, il est choisi par l'utilisateur et inscrit au journal.
    """

    notation: SequenceNotation
    library: ResidueLibrary | None = None
    name: str = "input.sequence"

    def run(self, data: str, context: StageContext) -> Peptide:
        peptide = parse_sequence(data, self.library, notation=self.notation)
        context.journal.record(
            "parse_sequence",
            "lecture de la séquence d'entrée",
            input=data,
            notation=peptide.metadata[SEQUENCE_NOTATION],
            n_terminus=peptide.n_terminus.key,
            c_terminus=peptide.c_terminus.key,
            length=len(peptide),
        )
        return peptide


@dataclass
class ProtonationStage:
    """Étape 2 : `Peptide` + pH de la configuration → `ProtonationResult`.

    Le pH vient de `PrepConfig.ph`, sans valeur par défaut. Les états imposés
    dans `PrepConfig.protonation` sont appliqués d'abord, puis chaque groupe
    ionisable reçoit son état (`pepforge.residues.protonation`). Chaque choix
    est journalisé avec son pKa et la source de ce pKa.
    """

    table: IonizableGroupTable | None = None
    name: str = "normalize.protonation"

    def run(self, data: Peptide, context: StageContext) -> ProtonationResult:
        config = context.config
        if config.ph is None:
            raise ConfigurationError(
                "pH de référence absent de la configuration (`ph`) : il est nécessaire pour "
                "choisir les états de protonation, et n'a pas de valeur par défaut"
            )
        table = self.table or load_ionizable_groups()
        peptide = self._impose(data, config.protonation, context)
        result = assign_protonation(peptide, config.ph, table)
        for assignment in result.assignments:
            context.journal.record(
                "assign_protonation",
                "état imposé par la configuration ou l'entrée"
                if assignment.imposed
                else f"forme majoritaire au pH (forme {assignment.form.label})",
                site=assignment.site,
                residue=assignment.residue,
                group=assignment.group,
                state=assignment.state,
                charge=assignment.charge,
                ph=result.ph,
                pka=assignment.pka,
                pka_source=table.sources[assignment.pka_source],
            )
        context.journal.record(
            "net_charge",
            "charge nette des états retenus et charge moyenne au pH",
            net_charge=result.net_charge,
            theoretical_net_charge=round(result.theoretical_net_charge, 4),
            ph=result.ph,
            data=result.data_source,
        )
        return result

    @staticmethod
    def _impose(peptide: Peptide, imposed: dict[str, str], context: StageContext) -> Peptide:
        for site, state in imposed.items():
            if not isinstance(state, str):
                raise ConfigurationError(f"`protonation[{site}]` : clé d'état attendue")
            if site == N_TERMINUS_SITE:
                peptide = peptide.with_n_terminus_protonation(state)
            elif site == C_TERMINUS_SITE:
                peptide = peptide.with_c_terminus_protonation(state)
            elif site.isdigit() and 1 <= int(site) <= len(peptide):
                peptide = peptide.with_protonation(int(site), state)
            else:
                raise ConfigurationError(
                    f"`protonation` : site {site!r} invalide (attendu « N », « C » ou une "
                    f"position de 1 à {len(peptide)})"
                )
            context.journal.record(
                "impose_protonation", "état imposé par la configuration", site=site, state=state
            )
        return peptide
