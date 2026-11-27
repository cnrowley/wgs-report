"""PharmCAT wrapper: star-allele calls and CPIC drug recommendations from the user's VCF.

PharmCAT (PharmGKB/CPIC, Apache-2.0) is not bundled. Download it yourself (see sources.PHARMCAT_RELEASE_URL) into a directory
that holds `pharmcat.jar` and the `pharmcat_pipeline` script (the "pipeline" tarball), and install its Python requirements
(`pip install colorama pandas packaging`). Java 17+ is needed.

Workflow: `pharmcat_pipeline -matcher` preprocesses the VCF (normalises it and extracts the PGx positions), then
`pharmcat.jar` runs the matcher, phenotyper and reporter on the preprocessed VCF. The jar is called directly because
only it accepts an outside call (`-po`) for genes PharmCAT cannot call from a short-read VCF, namely CYP2D6.

PharmCAT assumes positions absent from a variant-only VCF are reference. Read that as the same caveat as everywhere else in this
tool: confirm anything that matters (see the CRAM read-level checks).
"""
from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from .utils import log

TESTED_VERSIONS = ("3.4.0",)   # PharmCAT releases whose reporter JSON this parser was checked against (see tests/test_tools.py fixtures)

# Phenotypes that need no action on their own; drug rows are shown only for genes outside this set.
ROUTINE = re.compile(r"^(normal|no result|indeterminate|n/a|.*\bnormal (metabolizer|function)\b)", re.I)


def _all_reference(label: str) -> bool:
    """True for 'Reference/Reference' and 'rsN reference (C)/rsN reference (C)' style labels (nothing to act on)."""
    parts = [p.lower() for p in label.split("/")]
    return len(parts) == 2 and all("reference" in p and "variant" not in p for p in parts)


def _pharmcat_files(pharmcat_dir: Path) -> tuple[Path, Path]:
    d = Path(pharmcat_dir)
    jar, pipe = d / "pharmcat.jar", d / "pharmcat_pipeline"
    for p in (jar, pipe):
        if not p.exists():
            raise SystemExit(f"{p} not found. Unpack the PharmCAT 'pipeline' tarball into {d}.")
    return jar, pipe


def installed_version(jar: Path) -> str | None:
    """Ask the jar for its version (`--version`) before it is used; None if it cannot be read."""
    try:
        out = subprocess.run(["java", "-jar", str(jar), "--version"], capture_output=True, text=True, timeout=120, cwd=str(Path(jar).parent))
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(r"PharmCAT\s+v?(\d+\.\d+\.\d+\S*)", out.stdout + out.stderr)
    return m.group(1) if m else None


def check_install(jar: Path) -> str | None:
    """Fail early if Java is missing; warn (do not block) when the PharmCAT release is not one the parser was checked against."""
    if shutil.which("java") is None:
        raise SystemExit("PharmCAT needs Java 17 or newer on PATH (java -version). Install a JRE, or skip `wgsreport pharmcat`: "
                         "the built-in CPIC-style summary does not need it.")
    ver = installed_version(jar)
    if ver is None:
        log("PharmCAT: WARNING could not read the jar's version; output may not parse")
    elif ver not in TESTED_VERSIONS:
        log(f"PharmCAT: WARNING release {ver} is not one this parser was checked against ({', '.join(TESTED_VERSIONS)}); "
            "the report will carry a warning and the output is validated before use")
    return ver


def outside_call_file(path: Path, calls: dict[str, str]) -> Path:
    """Write PharmCAT's outside-call TSV (gene<TAB>diplotype), e.g. {'CYP2D6': '*1/*4'}."""
    path.write_text("".join(f"{g}\t{d}\n" for g, d in calls.items()))
    return path


def normalize_outside(calls: dict[str, str] | None) -> dict[str, str] | None:
    """PharmCAT does not know tandem notation with a *68 hybrid ('*68+x'). *68 is a non-functional CYP2D7 conversion, so the pair
    is passed as the allele that accompanies it."""
    if not calls:
        return calls
    return {g: (re.sub(r"\*68\+(?=\*)", "", d) if g == "CYP2D6" else d) for g, d in calls.items()}


def run(vcf: Path, out_dir: Path, pharmcat_dir: Path, outside_calls: dict[str, str] | None = None, sample: str = "sample",
        absent_to_ref: bool = True) -> Path:
    """Preprocess the VCF, run PharmCAT and return the path of the reporter JSON.

    absent_to_ref: treat PGx positions missing from a variant-only VCF as homozygous reference. PharmCAT marks this
    "dangerous" because a missing position can also mean no coverage; without it most genes return "No Result" and the
    matcher returns tied diplotypes. Confirm anything that matters with the CRAM read-level checks."""
    pharmcat_dir, vcf = Path(pharmcat_dir).resolve(), Path(vcf).resolve()
    jar, pipe = _pharmcat_files(pharmcat_dir)
    check_install(jar)
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base = re.sub(r"\.vcf(\.gz|\.bgz)?$", "", Path(vcf).name)
    log("PharmCAT: preprocessing the VCF (first run downloads its reference FASTA, which takes a few minutes)")
    # The pipeline script imports its sibling `pcat` package, so `python -I` cannot be used; run it from its own directory
    # with the interpreter of the tool's own virtualenv if there is one (third-party code stays out of our environment).
    py = next((q for q in (pharmcat_dir / "venv" / "bin" / "python", pharmcat_dir.parent / "venv" / "bin" / "python") if q.exists()), None)
    subprocess.run([str(py) if py else sys.executable, str(pipe), str(vcf), "-matcher", "-o", str(out_dir), "-bf", base] + (["--absent-to-ref"] if absent_to_ref else []),
                   check=True, cwd=str(jar.parent))
    pre = out_dir / f"{base}.preprocessed.vcf.bgz"
    cmd = ["java", "-jar", str(jar), "-vcf", str(pre), "-reporterJson", "-o", str(out_dir), "-bf", base]
    if outside_calls:
        cmd += ["-po", str(outside_call_file(out_dir / "outside_calls.tsv", normalize_outside(outside_calls)))]
    log("PharmCAT: matching alleles and building recommendations")
    subprocess.run(cmd, check=True)
    return out_dir / f"{base}.report.json"


def _label(d: dict) -> str:
    return d.get("label") or f"{d['allele1']['name']}/{d['allele2']['name']}" if d.get("allele1") and d.get("allele2") else d.get("label", "")


def summarize(report: dict, sources: tuple[str, ...] = ("CPIC Guideline Annotation",)) -> dict:
    """Reduce PharmCAT's reporter JSON to what the report shows: one row per gene, plus drug recommendations for
    genes whose phenotype is not routine. The full JSON stays on disk for anyone who wants every guideline."""
    if not isinstance(report.get("genes"), dict) or "drugs" not in report:
        raise ValueError("unexpected PharmCAT reporter JSON (no 'genes'/'drugs'); the output format has probably changed. "
                         f"This parser was checked against PharmCAT {', '.join(TESTED_VERSIONS)}. The built-in CPIC-style "
                         "summary (pgx.py) in the report does not depend on PharmCAT.")
    genes = []
    for gene, g in sorted((report.get("genes") or {}).items()):
        diplos = g.get("recommendationDiplotypes") or g.get("sourceDiplotypes") or []
        if not diplos:
            continue
        labels = [d.get("label") or _label(d) for d in diplos]
        phenos = sorted({p for d in diplos for p in (d.get("phenotypes") or [])}) or ["No result"]
        ambiguous = len(set(labels)) > 1
        diplotype = ("ambiguous: " + ", ".join(sorted(set(labels)))) if ambiguous else labels[0]
        genes.append({"gene": gene, "diplotype": diplotype, "phenotype": "; ".join(phenos), "ambiguous": ambiguous,
                      "activity_score": None if ambiguous else diplos[0].get("activityScore"), "call_source": g.get("callSource", ""),
                      "uncalled": list(g.get("uncalledHaplotypes") or []),
                      "undocumented_variation": bool(g.get("hasUndocumentedVariations")),
                      "actionable": any(not ROUTINE.search(p) for p in phenos) and not (not ambiguous and _all_reference(diplotype))})
    drugs, routine_n, seen = [], 0, set()
    by_gene = {x["gene"]: x for x in genes}
    for src in sources:
        for drug, rec in sorted((report.get("drugs") or {}).get(src, {}).items()):
            for guideline in rec.get("guidelines", []):
                for a in guideline.get("annotations", []):
                    gs = [list(k.keys())[0] for k in a.get("lookupKey", []) if isinstance(k, dict) and k]
                    text = html.unescape(a.get("drugRecommendation") or "")
                    if a.get("classification") == "No recommendation" or text.strip().lower() == "no recommendation":
                        continue
                    if not gs or not any(by_gene.get(x, {}).get("actionable") for x in gs):
                        routine_n += 1
                        continue
                    key = (drug, text)
                    if key in seen:
                        continue
                    seen.add(key)
                    drugs.append({"drug": drug, "genes": gs, "recommendation": text, "population": html.unescape(a.get("population") or ""),
                                  "implication": html.unescape("; ".join(a.get("implications") or [])),
                                  "strength": a.get("classification") or "", "source": src.split()[0]})
    ver = report.get("pharmcatVersion")
    warn = None if ver in TESTED_VERSIONS else (f"PharmCAT {ver} was not tested with this parser (tested: {', '.join(TESTED_VERSIONS)}); "
                                               "check the table against the full PharmCAT report.")
    return {"version": ver, "version_warning": warn, "data_version": report.get("dataVersion"),
            "genes": genes, "drugs": drugs, "routine_recommendations_hidden": routine_n}


def load_summary(path: Path) -> dict:
    return summarize(json.loads(Path(path).read_text()))
