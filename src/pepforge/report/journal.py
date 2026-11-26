"""Journal des opérations appliquées à une structure.

Toute opération qui modifie la structure enregistre : quoi, pourquoi, avec
quels paramètres.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from typing import Any

import yaml


@dataclass(frozen=True)
class JournalEntry:
    action: str
    reason: str
    parameters: dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(UTC).isoformat())


@dataclass
class Journal:
    entries: list[JournalEntry] = field(default_factory=list)

    def record(self, action: str, reason: str, **parameters: Any) -> JournalEntry:
        entry = JournalEntry(action=action, reason=reason, parameters=parameters)
        self.entries.append(entry)
        return entry

    def to_yaml(self) -> str:
        return yaml.safe_dump(
            [asdict(e) for e in self.entries], allow_unicode=True, sort_keys=False
        )
