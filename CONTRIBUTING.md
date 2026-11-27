# Contributing

Thank you for helping. Start with `CLAUDE.md` (rules and repo map) and `docs/EXTENDING.md` (how to add a module, a gene list or a score).

* Open an issue before a large change. Keep pull requests small and focused.
* Run `pytest` and `bash scripts/check_no_personal_data.sh` before every commit. **Never** include a real genome, report, `results.json`,
  sample ID, name or lab value in an issue, pull request or test.
* New predictions must meet the inclusion rule in `data/traits.yaml`; new thresholds must cite their source in the data file header.
* Use synthetic fixtures in tests. Network tests are marked `network` and run weekly in CI.
* Be honest in wording: each result states its limit, and nothing is presented as a diagnosis.
