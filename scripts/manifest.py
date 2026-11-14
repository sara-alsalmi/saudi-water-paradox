"""data/raw/manifest.csv — one row per downloaded file (idempotent upsert by local_path)."""
from __future__ import annotations

import csv
import hashlib
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

from common import ROOT

MANIFEST = ROOT / "data" / "raw" / "manifest.csv"
FIELDS = ["source", "local_path", "url", "version", "sha256", "size_bytes", "downloaded_utc", "notes"]
_LOCK = threading.Lock()   # downloads run in threads; serialise read-modify-write of the CSV


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read() -> dict[str, dict]:
    if not MANIFEST.exists():
        return {}
    with open(MANIFEST, newline="", encoding="utf-8") as f:
        return {r["local_path"]: r for r in csv.DictReader(f)}


def already_have(path: Path) -> bool:
    """True if file exists and matches the checksum recorded in the manifest."""
    rel = path.relative_to(ROOT).as_posix()
    with _LOCK:
        row = _read().get(rel)
    return bool(row) and path.exists() and sha256(path) == row["sha256"]


def record(source: str, path: Path, url: str, version: str, notes: str = "") -> None:
    rel = path.relative_to(ROOT).as_posix()
    new = {"source": source, "local_path": rel, "url": url, "version": version,
           "sha256": sha256(path), "size_bytes": path.stat().st_size,
           "downloaded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "notes": notes}
    with _LOCK:
        rows = _read()
        rows[rel] = new
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        tmp = MANIFEST.with_suffix(".csv.tmp")
        with open(tmp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for k in sorted(rows):
                w.writerow(rows[k])
        for attempt in range(10):          # Windows: target briefly locked while another reader has it open
            try:
                tmp.replace(MANIFEST)
                break
            except PermissionError:
                if attempt == 9:
                    raise
                time.sleep(0.5)
