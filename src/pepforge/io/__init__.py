"""Entrées et sorties : séquences, structures PDB et mmCIF, puis SMILES, MOL/SDF (à venir)."""

from pepforge.io.peptide_extraction import (
    ChainDecisions,
    Decision,
    ExtractionDecisions,
    ExtractionReport,
    ExtractionResult,
    analyze_chain,
    extract_chains,
    extract_peptides,
    structure_to_peptide,
)
from pepforge.io.retrieval import DownloadLimits, RetrievedFile, download_file, load_local_file
from pepforge.io.sequence import (
    SEQUENCE_NOTATION,
    SequenceNotation,
    convert_sequence,
    format_sequence,
    parse_sequence,
)
from pepforge.io.structure import (
    parse_structure,
    read_structure,
    read_structure_url,
    write_structure,
)
from pepforge.io.structure_formats import StructureFormat
from pepforge.io.structure_inventory import StructureInventory, inventory

__all__ = [
    "SEQUENCE_NOTATION",
    "ChainDecisions",
    "Decision",
    "DownloadLimits",
    "ExtractionDecisions",
    "ExtractionReport",
    "ExtractionResult",
    "RetrievedFile",
    "SequenceNotation",
    "StructureFormat",
    "StructureInventory",
    "analyze_chain",
    "convert_sequence",
    "download_file",
    "extract_chains",
    "extract_peptides",
    "format_sequence",
    "inventory",
    "load_local_file",
    "parse_sequence",
    "parse_structure",
    "read_structure",
    "read_structure_url",
    "structure_to_peptide",
    "write_structure",
]
