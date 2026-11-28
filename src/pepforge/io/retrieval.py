"""Obtention d'un fichier de structure : depuis le disque ou par téléchargement contrôlé.

Les deux modes d'importation demandés par Firas (27 novembre 2026) aboutissent
au même objet, `RetrievedFile` : le contenu exact reçu, son origine, sa taille,
son empreinte SHA-256 et sa date d'obtention. La lecture du format vient
ensuite (`pepforge.io.structure`) et ne sait pas d'où vient le fichier.

Téléchargement contrôlé (`download_file`) :

- schéma `https` seulement ; `http` n'est accepté que sur demande explicite
  (`DownloadLimits.allow_http`), et jamais après une redirection depuis `https` ;
- aucun identifiant dans l'URL (`https://nom:motdepasse@...` refusé) ;
- délai maximal pour chaque opération réseau (`timeout_s`) ;
- taille plafonnée **pendant** la lecture (`max_bytes`) : la longueur annoncée
  par le serveur est vérifiée d'abord, puis la lecture s'arrête au-delà du
  plafond, que l'annonce ait été honnête ou non ;
- réponse autre que 200 refusée ; réponse annoncée en HTML refusée (une page
  d'erreur ou de connexion n'est pas une structure) ;
- longueur reçue comparée à la longueur annoncée (téléchargement tronqué) ;
- fichier vide refusé.

Un fichier compressé par gzip (local ou téléchargé) est décompressé avec le
même plafond, appliqué à la taille décompressée. L'empreinte SHA-256 est celle
des octets **reçus**, avant décompression, pour pouvoir être comparée à celle
publiée par la source.

Rien n'est mis en cache et rien n'est écrit sur le disque : le contenu reste en
mémoire et n'est enregistré que si l'appelant le demande.
"""

from __future__ import annotations

import hashlib
import urllib.error
import urllib.parse
import urllib.request
import zlib
from dataclasses import dataclass
from datetime import UTC, datetime
from http.client import HTTPMessage
from pathlib import Path
from typing import IO

from pepforge.exceptions import RetrievalError

GZIP_MAGIC = b"\x1f\x8b"
_CHUNK = 1 << 16


@dataclass(frozen=True)
class DownloadLimits:
    """Bornes d'un téléchargement.

    Ce sont des réglages techniques, pas des choix scientifiques. Les valeurs
    par défaut couvrent une structure de récepteur ou de complexe (quelques
    mégaoctets) avec une large marge ; une entrée géante (ribosome, virus)
    demande de relever `max_bytes` explicitement.
    """

    timeout_s: float = 60.0
    """Délai maximal de chaque opération réseau (connexion, attente d'un bloc), en secondes."""

    max_bytes: int = 100 * 1024 * 1024
    """Taille maximale acceptée, reçue et décompressée, en octets (100 Mio)."""

    allow_http: bool = False
    """Accepter le schéma `http` non chiffré (serveur local de test, réseau interne)."""

    def __post_init__(self) -> None:
        if not self.timeout_s > 0:
            raise ValueError("timeout_s doit être strictement positif")
        if self.max_bytes <= 0:
            raise ValueError("max_bytes doit être strictement positif")


@dataclass(frozen=True)
class RetrievedFile:
    """Contenu d'un fichier obtenu, avec sa provenance."""

    content: bytes
    """Contenu prêt à lire (décompressé s'il était compressé)."""

    origin: str
    """Chemin absolu du fichier local, ou URL demandée."""

    mode: str
    """`local` ou `url`."""

    name: str
    """Nom du fichier (dernier élément du chemin ou de l'URL), qui porte l'extension."""

    sha256: str
    """Empreinte SHA-256 des octets reçus, avant décompression."""

    size: int
    """Nombre d'octets reçus, avant décompression."""

    retrieved_at: str
    """Date d'obtention, ISO 8601 en UTC."""

    compressed: bool = False
    final_url: str | None = None
    """URL après redirections (mode `url`)."""

    content_type: str | None = None
    """Type de contenu annoncé par le serveur (mode `url`)."""

    def provenance(self) -> dict[str, object]:
        """Ce que le journal enregistre de l'origine du fichier."""
        record: dict[str, object] = {
            "mode": self.mode,
            "origin": self.origin,
            "name": self.name,
            "sha256": self.sha256,
            "size_bytes": self.size,
            "retrieved_at": self.retrieved_at,
            "gzip": self.compressed,
        }
        if self.mode == "url":
            record["final_url"] = self.final_url
            record["content_type"] = self.content_type
        return record


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _decompress(raw: bytes, origin: str, max_bytes: int) -> bytes:
    decompressor = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    try:
        content = decompressor.decompress(raw, max_bytes + 1)
    except zlib.error as exc:
        raise RetrievalError(f"{origin} : contenu gzip illisible ({exc})") from None
    if len(content) > max_bytes or decompressor.unconsumed_tail:
        raise RetrievalError(
            f"{origin} : taille décompressée supérieure au plafond de {max_bytes} octets"
        )
    if not decompressor.eof:
        raise RetrievalError(f"{origin} : contenu gzip tronqué")
    if decompressor.unused_data:
        raise RetrievalError(
            f"{origin} : données après la fin du flux gzip (archive concaténée non prise en charge)"
        )
    return content


def _finish(
    raw: bytes, *, origin: str, mode: str, name: str, max_bytes: int, **extra: str | None
) -> RetrievedFile:
    if not raw:
        raise RetrievalError(f"{origin} : fichier vide")
    compressed = raw.startswith(GZIP_MAGIC)
    content = _decompress(raw, origin, max_bytes) if compressed else raw
    if not content.strip():
        raise RetrievalError(f"{origin} : fichier vide après décompression")
    return RetrievedFile(
        content=content,
        origin=origin,
        mode=mode,
        name=name,
        sha256=hashlib.sha256(raw).hexdigest(),
        size=len(raw),
        retrieved_at=_now(),
        compressed=compressed,
        final_url=extra.get("final_url"),
        content_type=extra.get("content_type"),
    )


def load_local_file(path: str | Path, max_bytes: int = DownloadLimits.max_bytes) -> RetrievedFile:
    """Lit un fichier local en entier ; absent, dossier ou vide : `RetrievalError`."""
    file_path = Path(path).expanduser()
    if not file_path.exists():
        raise RetrievalError(f"Fichier introuvable : {file_path}")
    if not file_path.is_file():
        raise RetrievalError(f"Ce n'est pas un fichier : {file_path}")
    size = file_path.stat().st_size
    if size > max_bytes:
        raise RetrievalError(
            f"{file_path} : {size} octets, au-delà du plafond de {max_bytes} octets"
        )
    raw = file_path.read_bytes()
    resolved = str(file_path.resolve())
    return _finish(raw, origin=resolved, mode="local", name=file_path.name, max_bytes=max_bytes)


class _GuardedRedirects(urllib.request.HTTPRedirectHandler):
    """Suit les redirections en refusant tout passage vers un schéma non autorisé."""

    def __init__(self, allowed: frozenset[str]) -> None:
        self.allowed = allowed

    def redirect_request(
        self,
        req: urllib.request.Request,
        fp: IO[bytes],
        code: int,
        msg: str,
        headers: HTTPMessage,
        newurl: str,
    ) -> urllib.request.Request | None:
        scheme = urllib.parse.urlsplit(newurl).scheme.lower()
        original = urllib.parse.urlsplit(req.full_url).scheme.lower()
        if scheme not in self.allowed or (original == "https" and scheme != "https"):
            raise RetrievalError(
                f"Redirection refusée de {req.full_url} vers {newurl} (schéma {scheme!r})"
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _check_url(url: str, limits: DownloadLimits) -> urllib.parse.SplitResult:
    parts = urllib.parse.urlsplit(url)
    allowed = {"https", "http"} if limits.allow_http else {"https"}
    scheme = parts.scheme.lower()
    if scheme not in allowed:
        hint = (
            " (http demande DownloadLimits(allow_http=True), --allow-http en ligne de commande)"
            if scheme == "http"
            else ""
        )
        raise RetrievalError(f"Schéma d'URL refusé : {scheme!r} dans {url!r}{hint}")
    if not parts.hostname:
        raise RetrievalError(f"URL sans nom d'hôte : {url!r}")
    if parts.username or parts.password:
        raise RetrievalError("URL refusée : elle contient un identifiant ou un mot de passe")
    return parts


def _name_from_url(parts: urllib.parse.SplitResult) -> str:
    return Path(urllib.parse.unquote(parts.path)).name


def download_file(url: str, limits: DownloadLimits | None = None) -> RetrievedFile:
    """Télécharge `url` sous les bornes de `limits` et renvoie le contenu vérifié."""
    limits = limits or DownloadLimits()
    parts = _check_url(url, limits)
    allowed = frozenset({"https", "http"} if limits.allow_http else {"https"})
    opener = urllib.request.build_opener(_GuardedRedirects(allowed))
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "pepforge-structure-import", "Accept-Encoding": "identity"},
    )
    try:
        with opener.open(request, timeout=limits.timeout_s) as response:
            status = response.status
            if status != 200:
                raise RetrievalError(f"{url} : réponse HTTP {status}, 200 attendu")
            content_type = response.headers.get("Content-Type")
            if content_type and content_type.split(";")[0].strip().lower() in (
                "text/html",
                "application/xhtml+xml",
            ):
                raise RetrievalError(
                    f"{url} : le serveur renvoie une page HTML ({content_type}), pas une structure"
                )
            encoding = (response.headers.get("Content-Encoding") or "identity").lower()
            if encoding not in ("identity", "gzip"):
                raise RetrievalError(
                    f"{url} : encodage de transfert non pris en charge {encoding!r}"
                )
            announced_text = response.headers.get("Content-Length")
            announced = int(announced_text) if announced_text and announced_text.isdigit() else None
            if announced is not None and announced > limits.max_bytes:
                raise RetrievalError(
                    f"{url} : {announced} octets annoncés, au-delà du plafond de "
                    f"{limits.max_bytes} octets"
                )
            raw = _read_capped(response, url, limits.max_bytes)
            final_url = response.geturl()
    except RetrievalError:
        raise
    except urllib.error.HTTPError as exc:
        raise RetrievalError(f"{url} : réponse HTTP {exc.code} ({exc.reason})") from None
    except urllib.error.URLError as exc:
        raise RetrievalError(f"{url} : téléchargement impossible ({exc.reason})") from None
    except TimeoutError:
        raise RetrievalError(f"{url} : délai de {limits.timeout_s} s dépassé") from None
    except OSError as exc:
        raise RetrievalError(f"{url} : erreur réseau ({exc})") from None
    if announced is not None and len(raw) != announced:
        raise RetrievalError(
            f"{url} : {len(raw)} octets reçus pour {announced} annoncés (téléchargement tronqué)"
        )
    return _finish(
        raw,
        origin=url,
        mode="url",
        name=_name_from_url(parts),
        max_bytes=limits.max_bytes,
        final_url=final_url,
        content_type=content_type,
    )


def _read_capped(response: IO[bytes], url: str, max_bytes: int) -> bytes:
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = response.read(_CHUNK)
        if not chunk:
            break
        total += len(chunk)
        if total > max_bytes:
            raise RetrievalError(
                f"{url} : plus de {max_bytes} octets reçus, téléchargement interrompu"
            )
        chunks.append(chunk)
    return b"".join(chunks)


__all__ = [
    "GZIP_MAGIC",
    "DownloadLimits",
    "RetrievalError",
    "RetrievedFile",
    "download_file",
    "load_local_file",
]
