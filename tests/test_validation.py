"""Failure-mode fixtures: each of these must break loudly if the scoring logic regresses (no real genomes involved)."""
import numpy as np
import pandas as pd
import pytest

from wgsreport import score
from wgsreport.pipeline import check_coverage
from .test_core import AF, pvar, sc


def test_allele_flip_inverts_weight_sign():
    """Same SNP reported with effect=ALT or effect=REF must give weights of opposite sign."""
    w_alt, _ = score.harmonize(sc([("1", 100, "G", "A", 0.5)]), pvar(), AF)
    w_ref, _ = score.harmonize(sc([("1", 100, "A", "G", 0.5)]), pvar(), AF)
    assert w_alt["1:100:A:G"] == 0.5 and w_ref["1:100:A:G"] == -0.5


def test_build_mismatch_is_detected():
    """Weights on shifted coordinates (e.g. a different build) match nothing, and prepare must refuse them."""
    shifted = sc([("1", 100 + 5000, "G", "A", 0.5), ("1", 200 + 5000, "T", "C", 0.2), ("2", 300 + 5000, "A", "C", 0.1)])
    w, st = score.harmonize(shifted, pvar(), AF)
    assert st["n_scored"] == 0 and st["coverage"] == 0
    with pytest.raises(SystemExit, match="genome build"):
        check_coverage("PGS000000", st)


def test_palindromic_ambiguous_snps_are_excluded():
    w, st = score.harmonize(sc([("2", 300, "T", "A", 1.0)]), pvar(), AF)   # A/T at MAF 0.5: strand cannot be resolved
    assert len(w[w != 0]) == 0 and st["n_palindromic_dropped"] == 1


def test_dosage_scaling_and_monotonic_extremes():
    p = pvar()
    w = pd.DataFrame({"s": [0.4, 0.3, 0.2, 0.1]}, index=p["id"])
    freq = np.full(4, 0.3)
    all_ref = score.user_scores(np.zeros(4), w, freq)["s"]
    het = score.user_scores(np.ones(4), w, freq)["s"]
    all_alt = score.user_scores(np.full(4, 2.0), w, freq)["s"]
    assert all_ref == 0 and het == pytest.approx(1.0) and all_alt == pytest.approx(2.0)   # 0/1/2 scale linearly
    assert all_ref < het < all_alt                                                         # extremes rank correctly
    neg = pd.DataFrame({"s": -w["s"]})
    assert score.user_scores(np.full(4, 2.0), neg, freq)["s"] < all_ref                    # protective weights reverse it


def test_every_snp_has_an_evidence_level():
    from wgsreport.registry import data_path
    ev = pd.read_csv(data_path("snp_evidence.tsv"), sep="\t", comment="#")
    panel = pd.read_csv(data_path("snp_panel.tsv"), sep="\t")
    assert set(panel["rsid"]) <= set(ev["rsid"]) and set(ev["level"]) <= {"A", "B", "C"}
