"""Compare measured lab values with the genetic tendency from the polygenic scores.

Input: a spreadsheet/CSV in "wide" form: first column = date, one column per test. A header row is the first row with
at least five text cells (a group-label row above it is ignored). Column names are matched case-insensitively against
ALIASES. Values must be in the units listed in TYPICAL (SI units typical of Canadian/European labs).

TYPICAL are *approximate* adult-male population means/SDs written from general knowledge, not from the user's lab
and not from the scores' training cohort; they are only used to put a measured value on a rough SD scale. Tests whose
values are analyser-dependent or strongly skewed are compared qualitatively only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

ALIASES = {
    "hemoglobin": "hb", "hemoglobin": "hb", "hgb": "hb", "hematocrit": "hct", "hematocrit": "hct", "rbc": "rbc",
    "mcv": "mcv", "mch": "mch", "rdw": "rdw", "wbc": "wbc", "platelets": "plt", "plt": "plt", "mpv": "mpv",
    "neutrophils": "neut", "lympocyctes": "lymph", "lymphocytes": "lymph", "monocytes": "mono", "hs-crp": "crp",
    "crp": "crp", "hdl": "hdl", "ldl": "ldl", "tg": "tg", "triglycerides": "tg", "lpa": "lpa", "apob": "apob",
}
# (mean, sd, unit) approximate adult-male; None = qualitative only
TYPICAL = {
    "hb": (150, 10, "g/L"), "hct": (0.44, 0.03, "fraction"), "rbc": (4.9, 0.4, "x10^12/L"), "mcv": (90, 4.5, "fL"),
    "mch": (30.5, 1.8, "pg"), "wbc": (6.3, 1.7, "x10^9/L"), "plt": (250, 55, "x10^9/L"), "neut": (3.6, 1.2, "x10^9/L"),
    "lymph": (1.9, 0.6, "x10^9/L"), "mono": (0.5, 0.15, "x10^9/L"), "hdl": (1.3, 0.3, "mmol/L"),
    "ldl": (3.2, 0.9, "mmol/L"), "apob": (0.95, 0.2, "g/L"),
    "rdw": None, "mpv": None, "tg": None, "crp": None, "lpa": None,
}

# Approximate typical adult reference intervals (low, high); a lab's own printed range takes precedence.
RANGES = {
    "hb": (130, 180), "hct": (0.40, 0.52), "rbc": (4.5, 5.9), "mcv": (80, 100), "mch": (27, 33), "wbc": (4.0, 11.0),
    "plt": (150, 400), "neut": (2.0, 7.5), "lymph": (1.0, 4.0), "mono": (0.2, 1.0),
}
# Canadian Cardiovascular Society 2021 dyslipidemia guideline (Pearson et al., CJC 2021;37:1129): in intermediate-risk adults a statin
# is recommended at LDL-C >= 3.5 mmol/L or ApoB >= 1.05 g/L; whether that applies depends on overall cardiovascular risk.
LIPID_THRESHOLDS = {"ldl": (3.5, "mmol/L", "LDL-C"), "apob": (1.05, "g/L", "ApoB")}

NOTES = {
    "lpa": ("Lipoprotein(a) is largely genetic and stable over life, and a measured value supersedes the polygenic estimate. The SNP "
            "score cannot see the KIV-2 repeat in LPA, the biggest single determinant, so a low score with a high measured value is "
            "possible. Check the unit on your own report: the Canadian Cardiovascular Society 2021 guideline lists >=100 nmol/L or "
            ">=50 mg/dL as separate thresholds, and there is no reliable universal conversion between them because particle mass "
            "depends on apo(a) isoform size. Read it alongside LDL/ApoB."),
    "ldl": ("If any of your LDL/ApoB results were measured on lipid-lowering treatment (see --treated-from), they reflect the "
            "treatment, not your genetic tendency; only untreated results are used in the comparison. Statins do not change Lp(a) much."),
    "apob": "ApoB measured on treatment reflects the treatment; compare the genetic tendency with untreated values.",
    "hdl": "HDL is strongly influenced by alcohol, exercise, smoking and diet, so measured values often differ from the genetic tendency.",
}


def load_wide(path) -> dict:
    p = str(path)
    raw = pd.read_csv(p, header=None) if p.lower().endswith(".csv") else pd.read_excel(p, header=None)
    hdr = next(i for i in range(len(raw)) if raw.iloc[i].map(lambda v: isinstance(v, str)).sum() >= 5)
    cols = [str(c).strip() if isinstance(c, str) else None for c in raw.iloc[hdr]]
    out: dict[str, list] = {}
    for i in range(hdr + 1, len(raw)):
        date = pd.to_datetime(raw.iloc[i, 0], errors="coerce")
        for j, c in enumerate(cols):
            if not c or j == 0:
                continue
            tid = ALIASES.get(c.lower())
            v = pd.to_numeric(raw.iloc[i, j], errors="coerce")
            if tid and pd.notna(v):
                out.setdefault(tid, []).append((None if pd.isna(date) else date.strftime("%Y-%m-%d"), float(v)))
    return out


TREATMENT_AFFECTED = {"ldl", "apob", "tg", "hdl", "lpa"}


def compare(traits: list[dict], measured: dict, treated_from: str | None = None) -> list[dict]:
    by_id = {t["id"]: t for t in traits if not t.get("skipped")}
    rows = []
    for tid, vals in measured.items():
        t = by_id.get(tid)
        if t is None or t["kind"] != "quantitative":
            continue
        typ = TYPICAL.get(tid)
        rng = RANGES.get(tid)
        in_range = bool(rng) and all(rng[0] <= v <= rng[1] for _, v in vals)
        untreated = [v for d, v in vals if not (treated_from and tid in TREATMENT_AFFECTED and tid != "lpa" and d and d >= treated_from)]
        mean_v = float(np.mean(untreated or [v for _, v in vals]))
        z_meas = None if not typ else (mean_v - typ[0]) / typ[1]
        exp = t["expected_sd"]
        if z_meas is None:
            verdict = "compare qualitatively"
        elif abs(exp) < 0.25:
            verdict = "genetics uninformative here"
        elif np.sign(exp) == np.sign(z_meas) or abs(z_meas) < 0.25:
            verdict = "directionally compatible"
        else:
            verdict = "differs from genetic tendency"
        if rng and in_range and verdict in ("directionally compatible", "differs from genetic tendency"):
            verdict += "; all values within the typical adult range, so this does not demonstrate a high or low value"
        flag = ""
        if tid in LIPID_THRESHOLDS and untreated:
            thr, unit, name = LIPID_THRESHOLDS[tid]
            if mean_v >= thr:
                flag = (f"Untreated {name} (mean {mean_v:.2f} {unit}) is at or above {thr:g} {unit}, the level at which the 2021 Canadian "
                        "guideline recommends a statin in intermediate-risk adults. Whether treatment applies depends on overall "
                        "cardiovascular risk, which this report does not assess.")
        if tid == "lpa" and untreated:
            hi = max(untreated)
            flag = (f"Reported as {hi:g} without a confirmed unit. If nmol/L, this is {hi / 100:.1f} times the 100 nmol/L threshold the Canadian "
                    "guideline lists as risk-enhancing. If mg/dL, it would be an exceptionally high value (the guideline's mg/dL "
                    "threshold is 50), so confirm the unit on the original report.") if hi >= 100 else ""
        rows.append({"id": tid, "label": t["label"], "values": vals, "clinical_flag": flag, "in_typical_range": in_range, "n_on_treatment": len(vals) - len(untreated), "mean": mean_v,
                     "unit": typ[2] if typ else "", "approx_z_measured": z_meas,
                     "genetic_percentile": t["percentile"], "expected_sd": exp, "verdict": verdict,
                     "note": NOTES.get(tid, "")})
    return rows
