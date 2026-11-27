"""One disclaimer, shown at the top of every output (reports, results.json, console)."""
TEXT = ("RESEARCH USE ONLY. This is not a medical or diagnostic test, was not produced in an accredited laboratory, and must not be used to make or change "
        "any medical decision. Polygenic scores describe statistical tendencies, not outcomes; they were developed mostly in European-ancestry cohorts and are "
        "less accurate, sometimes uninformative, in other ancestries (see the ancestry section). APOE genotype and a schizophrenia polygenic score are computed by "
        "default and can be sensitive: re-run with `--exclude apoe psychiatric` to leave them out.")


def html() -> str:
    t = TEXT.replace("`--exclude apoe psychiatric`", "<code>--exclude apoe psychiatric</code>")
    head, rest = t.split(". ", 1)
    return f"<div class='warn'><b>{head}.</b> {rest}</div>"
