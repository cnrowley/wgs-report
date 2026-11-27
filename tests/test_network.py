"""Upstream contract checks. Skipped by default; run with `pytest -m network` (also run weekly in CI)."""
import pytest

from wgsreport import doctor

pytestmark = pytest.mark.network


def test_upstream_endpoints_and_formats():
    bad = [r for r in doctor.check_network() if r[1] == "FAIL"]
    assert not bad, "\n".join(f"{n}: {m}" for n, _, m in bad)
