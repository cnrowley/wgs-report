"""`wgsreport doctor`: check that the tools, bundled data and upstream endpoints still look the way the code expects.

Run it after upgrading anything, when a download or parse fails, and from the weekly CI job. Each check returns
OK / WARN / FAIL with a message that says what to fix. Network checks are skipped with --offline.
"""
from __future__ import annotations

import gzip
import importlib
import re
import shutil
import subprocess
import urllib.request
import zlib

from . import sources
from .clinvar import acmg_genes
from .registry import load_registry, score_ids

TOOLS = {  # name: (version command, regex, minimum (major, minor), required?)
    "bcftools": (["bcftools", "--version"], r"bcftools (\d+)\.(\d+)", (1, 15), True),
    "samtools": (["samtools", "--version"], r"samtools (\d+)\.(\d+)", (1, 15), True),
    "plink2": (["plink2", "--version"], r"PLINK v(\d+)\.(\d+)", (2, 0), True),
    "bgzip": (["bgzip", "--version"], r"bgzip \(htslib\) (\d+)\.(\d+)", (1, 15), False),
    "tabix": (["tabix", "--version"], r"tabix \(htslib\) (\d+)\.(\d+)", (1, 15), False),
}
PY_PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "yaml"]


def _get(url: str, nbytes: int | None = None, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "wgsreport-doctor/0.1"})
    if nbytes:
        req.add_header("Range", f"bytes=0-{nbytes - 1}")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read(nbytes) if nbytes else r.read()


def check_tools():
    out = []
    for name, (cmd, rx, minimum, required) in TOOLS.items():
        if shutil.which(name) is None:
            out.append((f"tool:{name}", "FAIL" if required else "WARN", "not found on PATH"))
            continue
        try:
            txt = subprocess.run(cmd, capture_output=True, timeout=30)
            m = re.search(rx, (txt.stdout + txt.stderr).decode(errors="replace"))
            ver = (int(m.group(1)), int(m.group(2))) if m else None
        except Exception as e:  # noqa: BLE001
            out.append((f"tool:{name}", "WARN", f"could not read version: {e}"))
            continue
        if ver is None:
            out.append((f"tool:{name}", "WARN", "version output not recognised; the CLI may have changed"))
        elif ver < minimum:
            out.append((f"tool:{name}", "FAIL", f"{ver[0]}.{ver[1]} is older than required {minimum[0]}.{minimum[1]}"))
        else:
            out.append((f"tool:{name}", "OK", f"{ver[0]}.{ver[1]}"))
    for pkg in PY_PACKAGES:
        try:
            m = importlib.import_module(pkg)
            out.append((f"python:{pkg}", "OK", getattr(m, "__version__", "")))
        except ImportError:
            out.append((f"python:{pkg}", "FAIL", "not installed (pip install -e .)"))
    return out


def check_bundled_data():
    out = []
    reg = load_registry()
    ids = score_ids(reg)
    out.append(("data:traits", "OK", f"{len(reg['traits'])} traits, {len(ids)} scores"))
    n = len(acmg_genes())
    out.append(("data:acmg", "OK" if n >= 80 else "FAIL", f"{n} genes (v3.3 has 84)"))
    return out


def check_network():
    out = []

    def probe(name, fn):
        try:
            out.append((name, *fn()))
        except Exception as e:  # noqa: BLE001
            out.append((name, "FAIL", f"{type(e).__name__}: {e}"))

    def pgs():
        data = gzip.decompress(_get(sources.PGS_FTP.format(id="PGS000662")))
        head = [l for l in data.decode().splitlines() if not l.startswith("#")][0].split("\t")
        need = {"effect_allele", "effect_weight", "hm_chr", "hm_pos"}
        miss = need - set(head)
        return ("FAIL", f"harmonised file lacks columns {sorted(miss)}; update pgs.read_score") if miss else ("OK", "columns as expected")

    def clinvar():
        raw = _get(sources.CLINVAR_URL, nbytes=400_000)
        d = zlib.decompressobj(31).decompress(raw).decode(errors="ignore")
        ids = set(re.findall(r"##INFO=<ID=([A-Za-z0-9_]+)", d))
        miss = [f for f in sources.CLINVAR_INFO_FIELDS if f not in ids]
        return ("FAIL", f"INFO fields missing: {miss}; update clinvar.load_plp") if miss else ("OK", "INFO fields as expected")

    def ensembl():
        req = urllib.request.Request(sources.ENSEMBL_GFF_URL, method="HEAD", headers={"User-Agent": "wgsreport-doctor/0.1"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return ("OK", f"release {sources.ENSEMBL_RELEASE} reachable ({r.headers.get('Content-Length', '?')} bytes)")

    def ebi():
        seq = _get(sources.EBI_MD5_URL.format(md5=sources.CHRM_MD5)).decode().strip().replace("\n", "")
        return ("OK", "reference registry answers") if len(seq) == 16569 else ("FAIL", f"chrM length {len(seq)} != 16569")

    def pgs_rest():
        txt = _get(f"{sources.PGS_REST}/score/PGS000662?format=json").decode()
        return ("OK", "REST API answers") if '"id"' in txt else ("FAIL", "unexpected REST response")

    probe("net:pgs-catalog-file", pgs)
    probe("net:pgs-catalog-rest", pgs_rest)
    probe("net:clinvar", clinvar)
    probe("net:ensembl-gff", ensembl)
    probe("net:ebi-cram-registry", ebi)
    return out


def run_doctor(offline: bool = False) -> int:
    results = check_tools() + check_bundled_data() + ([] if offline else check_network())
    width = max(len(r[0]) for r in results)
    for name, status, msg in results:
        print(f"{status:5s} {name:<{width}}  {msg}")
    bad = [r for r in results if r[1] == "FAIL"]
    print(f"\n{len(bad)} failing check(s)" if bad else "\nall checks passed")
    return 1 if bad else 0
