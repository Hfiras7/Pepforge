"""Tests du journal des modifications."""

import yaml

from pepforge.report import Journal


def test_record_and_export() -> None:
    journal = Journal()
    journal.record("mutation", "test de l'effet de la charge", position=3, old="ASP", new="ASN")
    data = yaml.safe_load(journal.to_yaml())
    assert data[0]["action"] == "mutation"
    assert data[0]["parameters"]["new"] == "ASN"
    assert data[0]["timestamp"]
