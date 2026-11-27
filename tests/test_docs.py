"""Keeps the docs honest: every module is in the CLAUDE.md repo map, every CLI command is in the README, every data file has provenance."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "wgsreport"


def test_every_module_in_repo_map():
    claude = (ROOT / "CLAUDE.md").read_text()
    missing = [p.name for p in SRC.glob("*.py") if p.name != "__init__.py" and p.name not in claude and p.stem not in claude]
    assert not missing, f"add to the CLAUDE.md repo map: {missing}"


def test_every_command_in_readme():
    cli = (SRC / "cli.py").read_text()
    cmds = re.findall(r'sub\.add_parser\("([a-z0-9-]+)"', cli)
    readme = (ROOT / "README.md").read_text()
    missing = [c for c in cmds if c not in readme]
    assert cmds and not missing, f"document in README.md: {missing}"


GENERATED = {"snp_panel.tsv"}      # built by scripts/resolve_snp_panel.py from snp_panel_src.yaml, which carries the provenance


def test_data_files_have_provenance():
    for p in (SRC / "data").glob("*.tsv"):
        if p.name in GENERATED:
            continue
        head = p.read_text().splitlines()[:6]
        assert any(l.startswith("#") for l in head), f"{p.name} needs a '#' provenance header (source, date, how to rebuild)"
