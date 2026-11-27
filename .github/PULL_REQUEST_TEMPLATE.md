## What and why

## Checklist
- [ ] `pytest` passes; `bash scripts/check_no_personal_data.sh` is clean
- [ ] No personal data (genomes, reports, sample IDs, lab values, real tool output) in the diff or the description
- [ ] Thresholds, gene lists and category cut-offs come from a cited source recorded in the data file header (not from memory)
- [ ] Limits and "needs confirmation" wording are shown next to every new result; no absolute risks from percentiles
- [ ] New external tools are downloaded by the user (not bundled) and listed in `sources.py` and `docs/DATA_AND_LICENSES.md`
- [ ] Tests use synthetic fixtures and run offline
- [ ] README / `CLAUDE.md` repo map updated
