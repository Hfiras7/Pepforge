"""Tests de l'importation par URL : téléchargement contrôlé, sans Internet.

Un serveur HTTP local (127.0.0.1, port libre) sert les fichiers d'essai et
simule les réponses à refuser : erreur 404, page HTML, fichier trop gros,
longueur annoncée mensongère, téléchargement tronqué, serveur trop lent,
redirections. Le schéma `http` n'est accepté qu'avec `allow_http=True`, ce que
ces tests demandent explicitement.
"""

from __future__ import annotations

import gzip
import hashlib
import threading
import time
from collections.abc import Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from pepforge.exceptions import RetrievalError, StructureFormatError
from pepforge.io import DownloadLimits, download_file, read_structure, read_structure_url
from pepforge.io.retrieval import load_local_file
from pepforge.report.journal import Journal

DATA = Path(__file__).parent / "data" / "structures"
PDB_BYTES = (DATA / "capped_peptide.pdb").read_bytes()
CIF_BYTES = (DATA / "capped_peptide.cif").read_bytes()
GZ_BYTES = gzip.compress(PDB_BYTES, mtime=0)


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, *args: object) -> None:  # silence du serveur de test
        pass

    def _send(
        self, body: bytes, content_type: str = "text/plain", length: int | None = None
    ) -> None:
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body) if length is None else length))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        path = self.path
        if path == "/capped_peptide.pdb":
            self._send(PDB_BYTES)
        elif path == "/entry?id=x":
            self._send(CIF_BYTES, "chemical/x-mmcif")
        elif path == "/capped_peptide.pdb.gz":
            self._send(GZ_BYTES, "application/gzip")
        elif path == "/missing.pdb":
            self.send_error(404, "Not Found")
        elif path == "/login.pdb":
            self._send(b"<html><body>Connexion requise</body></html>", "text/html; charset=utf-8")
        elif path == "/big.pdb":
            self._send(PDB_BYTES * 10)
        elif path == "/lying.pdb":
            # Annonce une petite taille, envoie beaucoup plus, sans longueur fiable.
            self.send_response(200)
            self.send_header("Content-Type", "text/plain")
            self.end_headers()
            self.wfile.write(PDB_BYTES * 10)
            self.close_connection = True
        elif path == "/truncated.pdb":
            self._send(PDB_BYTES[:200], length=len(PDB_BYTES))
            self.close_connection = True
        elif path == "/slow.pdb":
            time.sleep(1.5)
            self._send(PDB_BYTES)
        elif path == "/partial.pdb":
            self.send_response(206)
            self.send_header("Content-Type", "text/plain")
            self.send_header("Content-Length", str(len(PDB_BYTES)))
            self.end_headers()
            self.wfile.write(PDB_BYTES)
        elif path == "/notes.pdb":
            self._send(b"Ceci n'est pas une structure.\n")
        elif path == "/empty.pdb":
            self._send(b"")
        elif path == "/redirect":
            self.send_response(302)
            self.send_header("Location", "/capped_peptide.pdb")
            self.end_headers()
        elif path == "/redirect-file":
            self.send_response(302)
            self.send_header("Location", "file:///etc/passwd")
            self.end_headers()
        else:
            self.send_error(404, "Not Found")


class _QuietServer(ThreadingHTTPServer):
    def handle_error(self, request: object, client_address: object) -> None:
        # Le client coupe volontairement certaines réponses (plafond, délai) : la
        # connexion rompue côté serveur est attendue.
        pass


@pytest.fixture(scope="module")
def server() -> Iterator[str]:
    httpd = _QuietServer(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{httpd.server_address[1]}"
    finally:
        httpd.shutdown()
        httpd.server_close()


LOCAL = DownloadLimits(allow_http=True, timeout_s=5)


def test_url_and_local_file_give_the_same_structure(server: str) -> None:
    journal = Journal()
    remote = read_structure_url(f"{server}/capped_peptide.pdb", journal=journal, limits=LOCAL)
    local = read_structure(DATA / "capped_peptide.pdb", journal=Journal())
    assert remote.atoms == local.atoms
    assert remote.connections == local.connections
    params = journal.entries[0].parameters
    assert params["mode"] == "url"
    assert params["origin"] == f"{server}/capped_peptide.pdb"
    assert params["sha256"] == hashlib.sha256(PDB_BYTES).hexdigest()
    assert params["size_bytes"] == len(PDB_BYTES)
    assert params["final_url"] == f"{server}/capped_peptide.pdb"


def test_url_without_extension_is_read_from_its_content(server: str) -> None:
    structure = read_structure_url(f"{server}/entry?id=x", journal=Journal(), limits=LOCAL)
    assert structure.metadata["structure_format"] == "mmcif"


def test_gzip_download_is_decompressed(server: str) -> None:
    retrieved = download_file(f"{server}/capped_peptide.pdb.gz", LOCAL)
    assert retrieved.compressed
    assert retrieved.content == PDB_BYTES
    # Empreinte des octets reçus (compressés), comparable à celle publiée par la source.
    assert retrieved.sha256 == hashlib.sha256(GZ_BYTES).hexdigest()
    assert retrieved.size == len(GZ_BYTES)


def test_redirect_is_followed_and_final_url_recorded(server: str) -> None:
    retrieved = download_file(f"{server}/redirect", LOCAL)
    assert retrieved.content == PDB_BYTES
    assert retrieved.final_url == f"{server}/capped_peptide.pdb"


def test_redirect_to_another_scheme_is_refused(server: str) -> None:
    # urllib refuse lui-même la redirection vers file:// ; l'erreur reste une RetrievalError.
    with pytest.raises(RetrievalError, match=r"Redirection refusée|not allowed"):
        download_file(f"{server}/redirect-file", LOCAL)


def test_redirect_from_https_to_http_is_refused_even_when_http_is_allowed() -> None:
    # Sans serveur TLS de test, la règle est vérifiée sur le gestionnaire de redirections.
    import urllib.request
    from email.message import Message

    from pepforge.io.retrieval import _GuardedRedirects

    handler = _GuardedRedirects(frozenset({"https", "http"}))
    request = urllib.request.Request("https://example.org/x.pdb")
    with pytest.raises(RetrievalError, match="Redirection refusée"):
        handler.redirect_request(
            request,
            None,
            302,
            "Found",
            Message(),
            "http://example.org/x.pdb",  # type: ignore[arg-type]
        )


def test_plain_http_is_refused_by_default(server: str) -> None:
    with pytest.raises(RetrievalError, match="allow_http"):
        download_file(f"{server}/capped_peptide.pdb")


@pytest.mark.parametrize(
    "url",
    ["ftp://example.org/x.pdb", "file:///etc/passwd", "x.pdb", "https:///x.pdb"],
)
def test_unsupported_urls_are_refused(url: str) -> None:
    with pytest.raises(RetrievalError):
        download_file(url, LOCAL)


def test_credentials_in_url_are_refused() -> None:
    with pytest.raises(RetrievalError, match="identifiant"):
        download_file("https://nom:secret@example.org/x.pdb")


def test_http_error_is_refused(server: str) -> None:
    with pytest.raises(RetrievalError, match="HTTP 404"):
        download_file(f"{server}/missing.pdb", LOCAL)


def test_partial_content_is_refused(server: str) -> None:
    with pytest.raises(RetrievalError, match="HTTP 206, 200 attendu"):
        download_file(f"{server}/partial.pdb", LOCAL)


def test_html_page_is_refused(server: str) -> None:
    with pytest.raises(RetrievalError, match="page HTML"):
        download_file(f"{server}/login.pdb", LOCAL)


def test_announced_size_above_the_cap_is_refused(server: str) -> None:
    limits = DownloadLimits(allow_http=True, timeout_s=5, max_bytes=len(PDB_BYTES) * 2)
    with pytest.raises(RetrievalError, match="annoncés, au-delà du plafond"):
        download_file(f"{server}/big.pdb", limits)


def test_body_above_the_cap_is_cut_even_without_announced_size(server: str) -> None:
    limits = DownloadLimits(allow_http=True, timeout_s=5, max_bytes=len(PDB_BYTES) * 2)
    with pytest.raises(RetrievalError, match="téléchargement interrompu"):
        download_file(f"{server}/lying.pdb", limits)


def test_truncated_download_is_refused(server: str) -> None:
    with pytest.raises(RetrievalError):
        download_file(f"{server}/truncated.pdb", LOCAL)


def test_timeout_is_enforced(server: str) -> None:
    limits = DownloadLimits(allow_http=True, timeout_s=0.3)
    with pytest.raises(RetrievalError, match=r"délai|téléchargement impossible|réseau"):
        download_file(f"{server}/slow.pdb", limits)


def test_empty_download_is_refused(server: str) -> None:
    with pytest.raises(RetrievalError, match="vide"):
        download_file(f"{server}/empty.pdb", LOCAL)


def test_downloaded_page_that_is_not_a_structure_is_refused(server: str) -> None:
    # Le serveur répond 200 en texte brut, mais le contenu n'est pas une structure.
    with pytest.raises(StructureFormatError, match="contenu non reconnu"):
        read_structure_url(f"{server}/notes.pdb", journal=Journal(), limits=LOCAL)


# --- Fichier local ---------------------------------------------------------------------


def test_local_file_above_the_cap_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "x.pdb"
    target.write_bytes(PDB_BYTES)
    with pytest.raises(RetrievalError, match="au-delà du plafond"):
        load_local_file(target, max_bytes=10)


def test_gzip_bomb_is_cut_at_the_cap(tmp_path: Path) -> None:
    target = tmp_path / "bomb.pdb.gz"
    target.write_bytes(gzip.compress(b"A" * 1_000_000))
    with pytest.raises(RetrievalError, match="décompressée"):
        load_local_file(target, max_bytes=100_000)


def test_truncated_gzip_is_refused(tmp_path: Path) -> None:
    target = tmp_path / "cut.pdb.gz"
    target.write_bytes(gzip.compress(PDB_BYTES)[:-12])
    with pytest.raises(RetrievalError, match="gzip"):
        load_local_file(target)


def test_directory_and_empty_file_are_refused(tmp_path: Path) -> None:
    with pytest.raises(RetrievalError, match="pas un fichier"):
        load_local_file(tmp_path)
    empty = tmp_path / "empty.pdb"
    empty.write_bytes(b"")
    with pytest.raises(RetrievalError, match="vide"):
        load_local_file(empty)


def test_limits_must_be_positive() -> None:
    with pytest.raises(ValueError):
        DownloadLimits(timeout_s=0)
    with pytest.raises(ValueError):
        DownloadLimits(max_bytes=0)
