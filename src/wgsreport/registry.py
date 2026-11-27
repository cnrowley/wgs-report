"""Load the trait registry shipped with the package."""
from __future__ import annotations

from importlib import resources
from pathlib import Path

import yaml


def data_path(name: str) -> Path:
    return Path(str(resources.files("wgsreport").joinpath("data", name)))


def load_registry(path: Path | None = None) -> dict:
    p = Path(path) if path else data_path("traits.yaml")
    with open(p) as f:
        reg = yaml.safe_load(f)
    seen = set()
    for t in reg["traits"]:
        if t["id"] in seen:
            raise ValueError(f"duplicate trait id {t['id']}")
        seen.add(t["id"])
        if t["kind"] == "quantitative" and "r2" not in t:
            raise ValueError(f"{t['id']}: quantitative traits need r2")
        if t["kind"] == "binary" and "or_per_sd" not in t:
            raise ValueError(f"{t['id']}: binary traits need or_per_sd")
        for g, alt in (t.get("alt_scores") or {}).items():
            if not {"pgs", "r2"} <= set(alt):
                raise ValueError(f"{t['id']}: alt_scores[{g}] needs pgs and r2")
    return reg


def score_ids(reg: dict) -> list[str]:
    ids = []
    for t in reg["traits"]:
        for k in ("pgs", "crosscheck", "variant_set"):
            if t.get(k) and t[k] not in ids:
                ids.append(t[k])
        for alt in (t.get("alt_scores") or {}).values():
            if alt["pgs"] not in ids:
                ids.append(alt["pgs"])
    return ids


def tier_from_r2(r2: float) -> str:
    """Same thresholds as the inclusion rule at the top of traits.yaml."""
    return "A" if r2 >= 0.20 else "B" if r2 >= 0.09 else "C"


def for_ancestry(t: dict, group: str, ancestry_ok: bool, available: set[str] | None = None) -> dict:
    """Pick the score and accuracy that fit the user's reference group (1000G super-population: EUR AFR AMR EAS SAS).

    Quantitative traits only. Evidence is per trait in traits.yaml: `r2_by_ancestry` (held-out R2 of the SAME score in other
    groups) and `alt_scores` (a different score that does better in a group, with its own R2). Rules, in order:
      1. not a clean single-group match, or EUR: the default score and its European figures (the report already warns for non-EUR);
      2. an `alt_scores` entry for the group whose R2 beats the default score's R2 in that group, and whose weights are in the panel;
      3. `r2_by_ancestry` for the group: same score, accuracy and tier recomputed from the group's R2 (tier C = percentile only);
      4. nothing published for the group: unchanged, flagged `portability` = 'unknown'.
    The returned dict is a copy with `portability` (text for the report) and `score_choice` set."""
    out = dict(t)
    out["score_choice"] = "default"
    if group == "EUR" or not ancestry_ok or t["kind"] != "quantitative":
        return out
    own = (t.get("r2_by_ancestry") or {}).get(group)
    alt = (t.get("alt_scores") or {}).get(group)
    if alt and (available is None or alt["pgs"] in available) and alt["r2"] > (own if own is not None else -1):
        out.update(pgs=alt["pgs"], r2=alt["r2"], r2_source=alt.get("source", ""), tier=tier_from_r2(alt["r2"]), score_choice="alt")
        out["portability"] = (f"Score switched to {alt['pgs']} because it was evaluated better in {group}-ancestry samples "
                              f"(R&sup2; {alt['r2']:.2f}" + (f" vs {own:.2f} for the default score" if own is not None else "") + ").")
    elif own is not None:
        out.update(r2=own, tier=tier_from_r2(own), r2_source=f"{t['pgs']} held-out in {group}-ancestry samples (see traits.yaml)", score_choice="same")
        out["portability"] = (f"Published accuracy of this score in {group}-ancestry samples is R&sup2; {own:.2f} "
                              f"(European {t['r2']:.2f}); tier recomputed.")
    else:
        out["portability"] = f"No published accuracy for this score in {group}-ancestry samples; the European figures overstate it."
    return out
