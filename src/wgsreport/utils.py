"""Small shared helpers: logging, external-tool handling, allele utilities."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def log(msg: str) -> None:
    print(f"[wgsreport {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def require_tools(*names: str) -> None:
    missing = [n for n in names if shutil.which(n) is None]
    if missing:
        raise SystemExit(
            "Missing required external tool(s): "
            + ", ".join(missing)
            + "\nInstall them (e.g. `conda env create -f environment.yml`) and retry."
        )


def run(cmd, **kwargs) -> subprocess.CompletedProcess:
    """Run a command (list or shell string) and raise with stderr on failure."""
    shell = isinstance(cmd, str)
    res = subprocess.run(cmd, shell=shell, text=True, capture_output=True,
                         executable="/bin/bash" if shell else None, **kwargs)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {cmd}\n{res.stderr[-3000:]}")
    return res


def complement(allele: str) -> str:
    return "".join(COMP[a] for a in allele)


def is_palindromic(a: str, b: str) -> bool:
    return len(a) == 1 and len(b) == 1 and COMP.get(a) == b


def download(url: str, dest: Path, retries: int = 5, timeout: int = 120) -> Path:
    """Download `url` to `dest` atomically with simple retry/backoff."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + f".{os.getpid()}.part")
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "wgsreport/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f)
            tmp.replace(dest)
            return dest
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"Failed to download {url}: {last}")
