# Extending wgs-report

This guide is for a person or an AI assistant adding a module, a data list or a score. Read `CLAUDE.md` first (hard rules, repo map, lessons
learned); this file is the step-by-step part.

## How the pieces fit
```
prepare (once)  ->  analyze (per person)  ->  report.html + results.json
                       ^ optional JSON summaries from separate commands
   pharmcat | hla-type | repeat-expansions | genome-qc | y-haplogroup | cram-check | copy-number
```
* Heavy or external work happens in its own command (`pharmcat`, `repeat-expansions`, ...) that writes a small **summary JSON**.
* `analyze` takes that JSON with a `--<name>-json` flag, stores it in `results.json` under one key, and `report.py` renders it with a
  `_<name>_html(summary)` function that returns `[]` when the summary is `None`. So a module never has to be run to produce a report.
* Everything that decides what a result *means* (thresholds, categories, wording) lives in a pure function that tests can call on a
  synthetic input. Wrappers that start programs stay thin.

## What to touch for what
| You want to add | Files |
|---|---|
| A tool-backed module (like repeat expansions) | `src/wgsreport/<name>.py`, `cli.py` (command + `analyze` flag), `analyze.py` (one line), `report.py` (`_<name>_html`), `sources.py` (URLs), `tests/test_<name>.py`, docs |
| A gene list or table | `src/wgsreport/data/<file>.tsv` with a provenance header, a `scripts/build_<x>.py` that regenerates it from the source, a loader in the module that uses it, a data test |
| A polygenic score | `data/traits.yaml` (see `CLAUDE.md` recipe), then `wgsreport prepare` |
| A curated SNP | `data/snp_panel_src.yaml`, run `scripts/resolve_snp_panel.py` |
| A report section only | `report.py` helper plus a test that renders it from a hand-built dict |

Fastest start: `python scripts/new_module.py <name>` writes a module skeleton and a passing test; then fill in the TODOs.

## Recipe: a tool-backed module
1. **Choose the source of truth.** Thresholds, gene lists and category cut-offs come from a cited file (STRchive, ClinGen, CPIC, a paper's
   table) loaded from `data/`, never typed from memory. Record the source and its date in the file header. A script rebuilds the file.
2. **Write `run(...)`** (starts the tool) and **`parse(...)`** (reads its output). Rules learned the hard way:
   * resolve every path with `Path(...).resolve()` before calling a tool that is run from another directory;
   * pass the non-interactive flag if the tool prompts (Yleaf `-force`), and turn off update checks and telemetry;
   * keep third-party code out of our environment (the tool's own venv or binary; `python -I` where the tool allows it);
   * delete personal intermediates (extracted reads, realigned BAMs) in a `finally` block.
3. **Write `summarize/classify(...)`** as a pure function: raw parsed rows in, `{rows, flagged, ...}` out. Say "not interpreted" when the
   input does not match what the thresholds describe (see `repeats.parse`, which only interprets a locus whose counted unit matches the
   pathogenic motif).
4. **Write `_<name>_html(...)`.** State the limits next to the result: what the method cannot see, which result needs confirmation, that a
   missing finding is not a negative clinical result. Never convert a score or percentile to an absolute risk.
5. **Wire it:** `cli.py` gets `<command>` (writes `<out>/<name>.json`) and `analyze --<name>-json`; `analyze.py` stores it in `out_json`;
   `report.py` calls the helper at the right place.
6. **Test with synthetic fixtures** shaped like the real output (`tests/test_tools.py` shows how). Cover the tie/ambiguous case, the
   "no result" case, and one false-lead the module could produce. Never commit real tool output.
7. **Run it for real on a public genome** (HG00096, or GIAB HG002 from 1000G/NIST) and compare with a second method if one exists.
8. **Docs:** one README paragraph with the command, the `CLAUDE.md` repo map line, `docs/DATA_AND_LICENSES.md` for any data or tool, and a
   `CLAUDE.md` "lessons" bullet for anything that surprised you.

## Definition of done (also in the PR template)
- [ ] `pytest` passes and `bash scripts/check_no_personal_data.sh` is clean; no personal data, sample IDs or real results in the diff.
- [ ] Every number a user sees has a stated source, threshold origin and limit. Ambiguity is shown, not resolved by taking the first option.
- [ ] Network is used only for public reference data and is listed in `sources.py`; tools are downloaded by the user, not bundled.
- [ ] Sensitive-but-valid results are computed by default with an `--exclude` switch, not hidden behind opt-in (see hard rules).
- [ ] The test for each new function runs offline in about a second.

## Working with Claude (or another assistant)
Start a session with: *"Read CLAUDE.md and docs/EXTENDING.md, then add <module> following the recipe. Use cited sources for every threshold,
write synthetic-fixture tests first, and run it on HG00096 before touching any personal data."* Useful habits:
* Ask it to print the real tool output once and build the parser from that, then reduce it to a synthetic fixture.
* Ask it to list false leads the new module can generate (for example a gene on a list for a different mechanism) and to test one.
* Keep personal VCF/CRAM outside the repository and give the assistant only the paths it needs; the repository's `.gitignore` and
  `scripts/check_no_personal_data.sh` are the safety net, not the plan.
* Set a local Git identity before committing (`git config user.name ...`; `git config user.email ...`) so that your real address does not
  end up in public history if you do not want it there.
