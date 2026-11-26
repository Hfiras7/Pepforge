"""Tests de l'encapsulation des outils externes, avec Python comme outil factice."""

import sys

import pytest

from pepforge.exceptions import ExternalToolError
from pepforge.tools import ExternalTool


def _python_tool(failure_patterns: tuple[str, ...] = ()) -> ExternalTool:
    return ExternalTool(executable=sys.executable, failure_patterns=failure_patterns)


def test_successful_run_records_version() -> None:
    result = _python_tool().run(["-c", "print('ok')"])
    assert result.stdout.strip() == "ok"
    assert result.version and "Python" in result.version


def test_nonzero_exit_raises() -> None:
    with pytest.raises(ExternalToolError):
        _python_tool().run(["-c", "import sys; sys.exit(2)"])


def test_failure_pattern_in_output_raises_even_with_zero_exit() -> None:
    tool = _python_tool(failure_patterns=(r"Could not find .*",))
    with pytest.raises(ExternalToolError, match="motif d'échec"):
        tool.run(["-c", "print('Could not find bond parameter for: CT - N')"])


def test_missing_executable() -> None:
    with pytest.raises(ExternalToolError, match="introuvable"):
        ExternalTool(executable="pepforge-outil-inexistant").path()
