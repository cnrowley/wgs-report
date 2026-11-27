"""Render the HTML report (self-contained: figures are embedded as base64 PNGs)."""
from __future__ import annotations

import base64
import html
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import norm  # noqa: E402

from . import disclaimer, hla_typing  # noqa: E402

E = html.escape
C_REF, C_USER, C_AXIS = "#8aa4c8", "#d1432b", "#555555"

CSS = """
:root{--bg:#fff;--fg:#1c2330;--mut:#5b6678;--line:#dfe4ec;--card:#f6f8fb;--warn:#fff4d6;--warnb:#e2b23a}
@media(prefers-color-scheme:dark){:root{--bg:#12161d;--fg:#e6eaf1;--mut:#9aa6b8;--line:#2a3140;--card:#1a2029;--warn:#3a3217;--warnb:#a17f1d}}
body{font:15px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--fg);max-width:1050px;margin:0 auto;padding:24px 16px}
h1{margin:.2em 0}h2{margin-top:2em;border-bottom:1px solid var(--line);padding-bottom:.2em}
table{border-collapse:collapse;width:100%;margin:.6em 0;font-size:14px}
th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top;overflow-wrap:break-word}
th{color:var(--mut);font-weight:600}.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 16px;margin:12px 0}
b.act-ok{color:#1a7f37}b.act-dose{color:#9a6700}b.act-bad{color:#cf222e}
.warn{background:var(--warn);border:1px solid var(--warnb);border-radius:8px;padding:10px 14px;margin:12px 0}
.mut{color:var(--mut)}.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(440px,1fr));gap:14px}
img{max-width:100%}.tag{display:inline-block;border:1px solid var(--line);border-radius:10px;padding:0 8px;font-size:12px;color:var(--mut)}
.big{font-size:22px;font-weight:650}@media(max-width:600px){.grid{grid-template-columns:1fr}}
"""


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _dist_fig(ref, user_raw, mu, sd, z, label):
    fig, ax = plt.subplots(figsize=(5.2, 2.0))
    zs = (np.asarray(ref) - mu) / sd
    ax.hist(zs, bins=30, color=C_REF, alpha=.85, density=True)
    x = np.linspace(-4, 4, 200)
    ax.plot(x, norm.pdf(x), color=C_AXIS, lw=1)
    ax.axvline((user_raw - mu) / sd, color=C_USER, lw=2)
    ax.set_xlim(-4.2, 4.2)
    ax.set_yticks([])
    ax.set_xlabel("polygenic score (SD vs reference population)", fontsize=8)
    ax.tick_params(labelsize=8)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title(label, fontsize=9, loc="left")
    return _png(fig)


def _pca_fig(panel_pcs, meta, user_pcs):
    df = panel_pcs.reset_index().merge(meta, on="IID")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    colors = {"AFR": "#e08a1e", "AMR": "#2f9e6e", "EAS": "#c7467e", "EUR": "#3b6fb6", "SAS": "#8a63c7"}
    for ax, (a, b) in zip(axs, (("PC1", "PC2"), ("PC3", "PC4"))):
        for sp, g in df.groupby("super_pop"):
            ax.scatter(g[a], g[b], s=5, alpha=.5, color=colors[sp], label=sp)
        ia, ib = int(a[2:]) - 1, int(b[2:]) - 1
        ax.scatter([user_pcs[ia]], [user_pcs[ib]], marker="*", s=240, color=C_USER, edgecolor="k", zorder=5)
        ax.set_xlabel(a, fontsize=9)
        ax.set_ylabel(b, fontsize=9)
        ax.tick_params(labelsize=8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", color=c, markersize=5, label=k) for k, c in colors.items()]
    handles.append(Line2D([], [], marker="*", ls="", color=C_USER, markeredgecolor="k", markersize=11, label="you"))
    axs[0].legend(handles=handles, fontsize=7)
    fig.tight_layout()
    return _png(fig)


def _ordinal(n: int) -> str:
    if 10 <= n % 100 <= 20:
        suf = "th"
    else:
        suf = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf}"


def _fmt_pct(p):
    if 1 <= p <= 99:
        return _ordinal(int(round(p)))
    return f"{p:.1f} (tail: imprecise)" if p < 100 else "above 99.9"


def _pct_range(z: float, n: int):
    """Approximate 95% sampling range for the percentile when the reference has n individuals (Gaussian model)."""
    from math import erf, sqrt
    se = sqrt(1.0 / n + z * z / (2.0 * n))
    cdf = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
    return 100 * cdf(z - 1.96 * se), 100 * cdf(z + 1.96 * se)


def _usable(n_scored: int, n_variants: int) -> str:
    f = n_scored / n_variants if n_variants else 0
    if f >= 1:
        return "100%"
    if f >= 0.9999:
        return ">99.99%"
    return f"{100 * f:.1f}%"


def _trait_card(r, fig_uri):
    if r.get("skipped"):
        return ""
    lines = [f"<div class='big'>{E(r['label'])} <span class='tag'>tier {r['tier']}</span></div>"]
    lo_p, hi_p = _pct_range(r["z"], r.get("n_ref", 500))
    shift = abs(r["z"] - r["z_raw"])
    lines.append(f"<div>Polygenic score: <b>{r['z']:+.2f} SD</b> &rarr; <b>{_fmt_pct(r['percentile'])} percentile</b> "
                 f"<span class='mut'>(plausible range {_ordinal(int(round(lo_p)))}&ndash;{_ordinal(int(round(hi_p)))} given the small reference; "
                 f"before ancestry adjustment: {r['z_raw']:+.2f} SD)</span></div>")
    if shift >= 0.5:
        lines.append(f"<div class='mut'>The ancestry adjustment moves this score by {shift:.2f} SD, so the result is sensitive to the "
                     "adjustment; read the percentile loosely.</div>")
    if r["kind"] == "quantitative" and r["tier"] == "C":
        lines.append(f"<div>Published accuracy of this score is low (R&sup2;&asymp;{r['r2']:.2f}); only the percentile is "
                     "reported, not a predicted trait value.</div>")
    elif r["kind"] == "quantitative":
        lines.append(f"<div>Typical value among people with a similar score: <b>{r['expected_sd']:+.2f} SD</b> from the population average. "
                     f"The score explains about {100 * r['r2']:.0f}% of the variation between people; most variation remains unexplained, "
                     "so this is not a prediction for an individual.</div>")
        if "abs_expected" in r:
            lines.append(f"<div>In units: typically about <b>{r['abs_expected']:.1f} {E(r['abs_unit'])}</b> for people with this score "
                         f"(reference mean {r['abs_mean']:g}, SD {r['abs_sd']:g}).</div>")
    elif r["tier"] == "C":
        lines.append("<div>Published accuracy of this score is low; only the percentile is reported.</div>")
    else:
        lines.append(f"<div>Odds relative to an average-score person: <b>&times;{r['relative_odds']:.2f}</b> "
                     f"<span class='mut'>(OR per SD = {r['or_per_sd']}). This is a relative odds multiplier, not a change in absolute risk "
                     "or in percent probability.</span></div>")
    if "crosscheck_z" in r:
        lines.append(f"<div class='mut'>Independent cross-check score {r['crosscheck_pgs']}: {r['crosscheck_z']:+.2f} SD</div>")
    if r.get("portability"):
        lines.append(f"<div class='warn'>{r['portability']}</div>")
    if "variant_set" in r:
        lines.append(_variant_set_html(r["variant_set"]))
    lines.append(f"<div class='mut'>{r['n_scored']:,}/{r['n_variants']:,} score variants usable "
                 f"({_usable(r['n_scored'], r['n_variants'])}) &middot; {E(r['pgs'])} &middot; {E(r['cite'])}</div>")
    if r.get("low_coverage"):
        lines.append("<div class='warn'>Low variant coverage; treat this score as unreliable.</div>")
    if r.get("validation"):
        lines.append("<details><summary>Published validation (PGC3 supplement)</summary><ul>" + "".join(f"<li>{E(v)}</li>" for v in r["validation"]) + "</ul></details>")
    if fig_uri:
        lines.append(f"<img src='{fig_uri}' alt='score distribution'>")
    return "<div class='card'>" + "\n".join(lines) + "</div>"


def _short(text, n=3, sep="|"):
    """Keep the first n entries of a long '|'-separated ClinVar field."""
    parts = [p.strip() for p in str(text).replace("_", " ").split(sep) if p.strip()]
    return "; ".join(parts[:n]) + (f" (+{len(parts) - n} more)" if len(parts) > n else "")



ACTION_CLASS = {"use normally": "act-ok", "adjust dose": "act-dose", "avoid / choose alternative": "act-bad", "may be less effective": "act-dose", "confirm first": "act-bad"}


def _pgx_consequence_html(rows) -> str:
    """What each result means for each drug: faster/slower metabolism, more/less effect, and use normally / adjust / avoid."""
    body = []
    for x in rows:
        for c in x.get("consequences") or []:
            body.append(f"<tr><td>{E(x['gene'])}</td><td>{E(c['drug'])}</td><td>{E(c['effect'])}</td>"
                        f"<td><b class='{ACTION_CLASS.get(c['action'], '')}'>{E(c['action'])}</b></td></tr>")
    if not body:
        return ""
    return ("<h3>What each result means for your medicines</h3><table><tr><th>Gene</th><th>Drug</th><th>Consequence for you</th>"
            "<th>Usual advice</th></tr>" + "".join(body) + "</table><p class='mut'>Summarised from CPIC guidelines for the alleles tested; "
            "this is general information, not dosing advice. &ldquo;Use normally&rdquo; means no tested variant changes the usual advice, not that every variant was tested.</p>")


def _pharmcat_html(pc) -> list:
    if not pc:
        return []
    out = [f"<h3>PharmCAT star-allele calls and CPIC recommendations (PharmCAT {E(str(pc.get('version') or ''))})</h3>"
           "<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Activity score</th><th>Call from</th></tr>"]
    if pc.get("version_warning"):
        out.insert(1, f"<div class='warn'>{E(pc['version_warning'])}</div>")
    no_result = []
    for g in pc.get("genes", []):
        if g["diplotype"].startswith("Unknown") or g["phenotype"].lower().startswith("no result") and g["diplotype"].lower().startswith("unknown"):
            no_result.append(g["gene"])
            continue
        flag = " <b>&bull;</b>" if g.get("actionable") else ""
        note = "outside caller" if g.get("call_source") == "OUTSIDE" else "VCF"
        out.append(f"<tr><td>{E(g['gene'])}{flag}</td><td>{E(g['diplotype'])}</td><td>{E(g['phenotype'])}</td>"
                   f"<td>{E(str(g['activity_score'])) if g.get('activity_score') is not None else ''}</td><td>{E(note)}</td></tr>")
    out.append("</table>")
    drugs = pc.get("drugs", [])
    if drugs:
        out.append("<p><b>&bull;</b> marks genes with a non-routine phenotype. CPIC recommendations for those genes:</p>"
                   "<table><tr><th>Drug</th><th>Gene</th><th>Recommendation (CPIC)</th><th>Clinical setting</th><th>Strength</th></tr>" +
                   "".join(f"<tr><td>{E(d['drug'])}</td><td>{E(', '.join(d['genes']))}</td><td>{E(d['recommendation'])}</td>"
                           f"<td>{E(d.get('population') or '')}</td><td>{E(d['strength'])}</td></tr>" for d in drugs) + "</table>")
    if no_result:
        out.append(f"<p class='mut'>No result from the VCF for: {E(', '.join(no_result))} (HLA-A/B are typed from reads in the HLA section below).</p>")
    out.append(f"<p class='mut'>Calls from the VCF assume that PGx positions absent from a variant-only VCF are reference. {pc.get('routine_recommendations_hidden', 0)} CPIC recommendations for routine phenotypes are not shown (they say "
               "to use standard dosing). Star-allele calls from a variant-only VCF cannot see structural variants or untested alleles, "
               "so &ldquo;normal&rdquo; means <i>no tested non-reference allele found</i>. These are research results: "
               "do not change any medicine on the strength of this page. Full JSON from PharmCAT stays on disk.</p>")
    return out


def _hla_html(h) -> list:
    if not h:
        return []
    out = ["<h3>HLA typing from reads (T1K)</h3><table><tr><th>Gene</th><th>Alleles</th><th>Zygosity</th><th>Quality</th></tr>"]
    for c in [x for x in h.get("calls", []) if x["gene"] in hla_typing.CLASSICAL]:
        extra = " (tied with other alleles)" if c.get("tied") else ""
        q = "" if c.get("quality") is None else f"{c['quality']:g}" + (" (low)" if c.get("low_confidence") else "")
        out.append(f"<tr><td>{E(c['gene'])}</td><td>{E(', '.join(c['alleles']) or 'not typed')}{E(extra)}</td><td>{E(c['zygosity'])}</td><td>{E(q)}</td></tr>")
    out.append("</table>")
    fl = h.get("flagged", [])
    if fl:
        out.append("<p><b>Alleles with guideline-level clinical associations:</b></p><table><tr><th>Allele</th><th>Copies</th><th>Meaning</th><th>Source</th></tr>" +
                   "".join(f"<tr><td>{E(f['allele'])}</td><td>{f['copies']}</td><td>{E(f['meaning'])}</td><td>{E(f['source'])}</td></tr>" for f in fl) + "</table>")
    else:
        out.append("<p>None of B*57:01, B*15:02, A*31:01, B*58:01 or B*27:05 detected.</p>")
    if h.get("celiac_typed"):
        c = h["celiac_compatible"]
        hits = [k for k, v in c.items() if v]
        out.append("<p><b>Celiac-associated DQ haplotypes:</b> " +
                   (f"alleles compatible with {E(', '.join(hits))}" if hits else "DQ2.5, DQ2.2 and DQ8 allele pairs not present") +
                   ". DQA1 and DQB1 phase (cis vs trans) cannot be read from short reads, so a compatible pair may not sit on one haplotype. "
                   "DQ2.2 alone carries little risk (the risk is mainly DQ2.5 and DQ8). Most carriers never develop celiac disease.</p>")
    out.append("<p class='mut'>Typing from 30x short reads is usually right at two-field resolution for the common classical genes and weaker for rare "
               "alleles and for DQ/DP/DRB genes. Confirm any result that matters (above all a drug hypersensitivity allele) with a clinical HLA test.</p>")
    return out


def _qc_html(q) -> list:
    if not q:
        return []
    out = ["<h2>Sex chromosomes, chromosome dosage and relatedness of your parents</h2>"]
    d = q.get("dosage")
    if d:
        r = d["ratios"]
        out.append(f"<p>Sex-chromosome dosage: <b>{E(d['sex_chromosome_call'])}</b> (chrX depth {r['chrX']:.2f}&times;, chrY {r['chrY']:.2f}&times; "
                   f"the autosomal median). {E(d['note'])}</p>")
        out.append("<p>" + ("Whole-chromosome gains or losses: <b>none</b> (every autosome within 0.65&ndash;1.25&times; the median; "
                            f"range {min(v for c, v in r.items() if c not in ('chrX', 'chrY')):.2f}&ndash;{max(v for c, v in r.items() if c not in ('chrX', 'chrY')):.2f})."
                            if not d["autosome_flags"] else "Possible whole-chromosome gain or loss: <b>" + E("; ".join(d["autosome_flags"])) +
                            "</b>. Confirm with a clinical karyotype or microarray.") + "</p>")
    xh = q.get("x_het")
    if xh and xh.get("het_fraction") is not None:
        out.append(f"<p>Heterozygous fraction of non-pseudoautosomal chrX calls: {xh['het_fraction']:.1%} ({xh['n_calls']:,} calls). "
                   "Near zero fits one X chromosome; about a third fits two.</p>")
    r = q.get("roh")
    if r:
        f = r["f_roh"]
        out.append(f"<p>Runs of homozygosity of at least {r['min_run_bp'] / 1e6:g} Mb: <b>{r['n_runs']}</b>, totalling {r['total_bp'] / 1e6:.1f} Mb "
                   f"(<b>{f:.2%}</b> of the autosomal genome; longest {r['longest_bp'] / 1e6:.1f} Mb). For comparison, the expected inbreeding coefficient for "
                   + ", ".join(f"{E(c['relationship'])} {c['expected_f']:.2%}" for c in r["context"]) +
                   " (theory; real values vary). Long runs mean both copies of that stretch came from a shared ancestor, which matters for recessive conditions. "
                   "Most people with unrelated parents have well under 1&ndash;2%. Long runs can also come from small ancestral populations.</p>")
    out.append("<p class='mut'>Dosage uses sampled read depth and detects whole-chromosome changes only (not mosaicism below roughly 20%, small copy-number "
               "changes or balanced rearrangements). Runs of homozygosity are inferred from genotypes without allele frequencies.</p>")
    return out


def variant_reading(r: dict, related_carried: list[str]) -> tuple[str, str] | None:
    """Headline for a trait that a few large-effect variants decide (Lp(a)): (reading, value) or None to use the polygenic score.

    Used when the trait sets `variant_set_primary`. The explicit-variant total (units of the trait) or a carried `related_rsids`
    variant outranks the genome-wide percentile, which describes only the diffuse background."""
    v = r.get("variant_set")
    if not v or not r.get("variant_set_primary"):
        return None
    unit, thr, tot = E(v.get("unit", "")), r.get("variant_set_threshold"), v["total"]
    name = E(r["label"].split(" (")[0])
    if related_carried or (thr is not None and tot >= thr):
        why = []
        if related_carried:
            why.append("carries " + ", ".join(E(x) for x in related_carried))
        if tot > 0:
            why.append(f"explicit-variant score adds {tot:+.0f} {unit}")
        return (f"<b>genetic variants predict raised {name}</b> ({'; '.join(why)}; a contribution above a non-carrier, not a total concentration)",
                f"&ge; {tot:.0f} {unit} above a non-carrier" if tot > 0 else "raised")
    if tot > 0:
        return (f"modest genetic contribution to {name} ({tot:+.0f} {unit} from carried variants; not a total concentration)", f"+{tot:.0f} {unit}")
    return None


def _variant_set_html(v: dict) -> str:
    unit = E(v.get("unit", ""))
    top = "; ".join(f"{E(x.get('rsid') or '%s:%s' % (x['chrom'], x['pos']))} {E(x['effect_allele'])} &times;{x['copies']} ({x['contribution']:+.1f})"
                    for x in v["carried"][:6]) or "none of the variants"
    ref = v.get("reference")
    rank = (f" Among {ref['n']} {E(ref['group'])} reference genomes, {100 * ref['fraction_at_or_above']:.0f}% have a total this high or higher "
            f"and {100 * ref['fraction_any_carried']:.0f}% carry any of these variants." if ref else "")
    return (f"<div>Explicit-variant cross-check {E(v['pgs'])} ({v['n_variants']} variants, weights in {unit}): carried variants add "
            f"<b>{v['total']:+.0f} {unit}</b> to a non-carrier. Carried: {top}.{rank} "
            f"<span class='mut'>Additive per copy of the effect allele, so a heterozygote counts once. Not standardised against a reference panel; "
            f"{v['assumed_reference']} sites absent from the VCF are assumed reference.</span></div>")


def _repeat_note(gene: str, category: str) -> str:
    if category == "pathogenic range":
        return ("Length in the disease-associated range. ExpansionHunter's size is an estimate and only a lower bound beyond the read length: "
                "confirm with a clinical repeat-primed PCR or long-read assay before drawing any conclusion.")
    note = "Intermediate (gray-zone) length: usually not expected to cause the disease in the carrier. "
    if gene == "HTT":
        return note + ("Alleles of this size can be unstable and expand into the disease range when passed on, mainly through paternal "
                       "transmission; the chance per child is low and rises with allele length.")
    return note + ("Such alleles can be unstable when passed on. Interruptions in the repeat, which ExpansionHunter does not see, often make "
                   "alleles of this length harmless, so a pure-repeat count overstates the concern.")


def _repeats_html(rp) -> list:
    if not rp:
        return []
    out = ["<h2>Repeat expansions (ExpansionHunter)</h2>"]
    fl = rp.get("flagged", [])
    interp = [x for x in rp.get("loci", []) if x.get("interpreted")]
    out.append(f"<p>{len(interp)} disease loci were genotyped and compared with the pathogenic-size thresholds from STRchive; "
               f"<b>{len(interp) - len(fl)}</b> are in the normal range"
               + (f" and <b>{len(fl)}</b> {'is' if len(fl) == 1 else 'are'} outside it." if fl else ", none outside it.") + "</p>")
    if fl:
        out.append("<table><tr><th>Gene</th><th>Disease</th><th>Repeat unit</th><th>Your repeat counts</th><th>Category</th><th>Normal up to</th><th>Meaning</th></tr>")
        for x in fl:
            out.append(f"<tr><td>{E(x['gene'])}</td><td>{E(x['disease'])} ({E(x['inheritance'])})</td><td>{E(x['unit'])}</td>"
                       f"<td>{E(x['genotype'])}{' (CI ' + E(x['ci']) + ')' if len(set(x['ci'].replace('/', '-').split('-'))) > 2 else ''}{' &ndash; lower bound' if x.get('beyond_read_length') else ''}</td>"
                       f"<td><b>{E(x['worst'])}</b></td><td>{x['thresholds']['benign_max']}</td><td>{E(_repeat_note(x['gene'], x['worst']))}</td></tr>")
        out.append("</table>")
    ok = ", ".join(sorted({x["gene"] for x in interp if x["worst"] == "normal range"}))
    if ok:
        out.append(f"<p class='mut'>Normal range: {E(ok)}.</p>")
    out.append("<p class='mut'>Screen of about 30 known loci only; sizes are reliable up to about the read length, and a large expansion appears as a lower "
               "bound. Reads whose mate maps to an alternate contig absent from the reference are not paired, and repeat interruptions are not "
               "assessed. A result outside the normal range is a lead for a clinical test, not a diagnosis.</p>")
    return out

def _table(df, cols, heads=None):
    if df is None or len(df) == 0:
        return "<p class='mut'>None.</p>"
    heads = heads or cols
    out = ["<table><tr>" + "".join(f"<th>{E(h)}</th>" for h in heads) + "</tr>"]
    for _, row in df.iterrows():
        out.append("<tr>" + "".join(f"<td>{E(str(row[c]))}</td>" for c in cols) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


CANCER_PGS = ("prostate", "melanoma")


def _cancer_html(res: dict, cv_hits) -> list:
    """Inherited cancer susceptibility: monogenic findings and polygenic scores kept apart, with what was and was not assessed."""
    import pandas as pd
    parts = ["<h2>Inherited cancer susceptibility</h2>",
             "<p>Two different things are shown separately and are <b>not</b> combined into one risk: rare pathogenic variants in "
             "predisposition genes (large effects, few people) and polygenic scores (many small effects, everyone has a value). "
             "No absolute lifetime risk is calculated: that needs family history, age, sex and a validated model.</p>"]
    has = cv_hits is not None and len(cv_hits) and "cancer_gene" in cv_hits
    hits = cv_hits[cv_hits["cancer_gene"].astype(str) != ""] if has else pd.DataFrame()
    if len(hits):
        hits = hits.assign(af_1000g=hits["af_1000g"].map(lambda x: "" if pd.isna(x) else f"{x:.4f}"))
        dominant = hits["cancer_gene"].astype(str).str.contains(r"AD|XL|SD", regex=True) | (hits["pattern"] != "single heterozygous")
        recessive_carrier = hits[~dominant]
        hits = hits[dominant]
    else:
        recessive_carrier = hits
    parts.append("<h3>A. Pathogenic / likely pathogenic variants in cancer-predisposition genes</h3>")
    if len(hits):
        t = hits.copy()
        t["condition"] = t["condition"].map(lambda x: _short(x, 3))
        t["gene (inheritance)"] = t["gene"] + " (" + t["cancer_gene"].astype(str) + ")"
        t["significance"] = t["significance"].map(lambda x: _short(x, 2, "/"))
        parts.append(_table(t, ["gene (inheritance)", "rsid", "zygosity", "pattern", "condition", "significance", "stars", "consequence", "af_1000g"]))
        parts.append("<p class='mut'>Inheritance is the mode curated for the gene (a gene can be dominant for cancer risk and recessive for another disease; "
                     "one copy in a recessive-only gene is carrier status). Penetrance depends on the variant, sex, family history and ancestry and is "
                     "not shown here because no number is embedded without a verified citation: use GeneReviews/ClinGen for the gene.</p>")
    else:
        parts.append("<p>None found. <b>Not</b> a negative clinical result (see C).</p>")
    if len(recessive_carrier):
        rc = recessive_carrier.copy()
        rc["condition"] = rc["condition"].map(lambda x: _short(x, 3))
        rc["gene (inheritance)"] = rc["gene"] + " (" + rc["cancer_gene"].astype(str) + ")"
        parts.append("<h4>Single copy in a recessive-only cancer-predisposition gene</h4><p class='mut'>Carrier status: one copy is not expected to raise "
                     "your own cancer risk (some genes, for example MUTYH or BLM, have weaker or debated heterozygote effects); relevant for family planning.</p>"
                     + _table(rc, ["gene (inheritance)", "rsid", "zygosity", "condition", "significance", "stars", "consequence", "af_1000g"]))
    low = pd.DataFrame(res.get("cancer_low_review") or [])
    if len(low):
        low = low.copy()
        low["af_1000g"] = low["af_1000g"].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
        low["condition"] = low["condition"].map(lambda x: _short(x, 3))
        low["gene (inheritance)"] = low["gene"] + " (" + low["cancer_gene"].astype(str) + ")"
        parts.append("<h4>Lower review status (ClinVar P/LP with fewer than 2 review stars)</h4>"
                     "<p class='mut'>Single-submitter or unreviewed classifications: weaker evidence, shown so they are not silently missed. "
                     "Treat as leads to check against ClinVar, not as findings.</p>"
                     + _table(low, ["gene (inheritance)", "rsid", "zygosity", "pattern", "condition", "significance", "stars", "consequence", "af_1000g"]))
    parts.append("<h3>B. Polygenic scores</h3>")
    rows = [r for r in res.get("traits", []) if r.get("id") in CANCER_PGS]
    if rows:
        li = "".join(f"<li><b>{E(r['label'])}</b>: tier {E(str(r['tier']))}, "
                     + ("not computed (sex-specific score)" if r.get("skipped") else
                        f"percentile {r.get('percentile', float('nan')):.0f}" if r.get("percentile") is not None else "see the disease-risk table")
                     + f"; reference {E(str(res.get('reference_group', '')))}.</li>" for r in rows)
        parts.append(f"<ul>{li}</ul><p class='mut'>Percentile in a reference of about 500 people, not a probability of cancer; accuracy is the published "
                     "held-out value for the ancestry it was tested in, and is lower elsewhere. Full details are in the disease-risk section.</p>")
    else:
        parts.append("<p>No cancer polygenic score was computed for this genome (prostate and melanoma are the ones that meet the inclusion rule).</p>")
    parts.append("<h3>C. What was and was not assessed</h3><ul>"
                 "<li><b>Assessed</b>: ClinVar pathogenic / likely pathogenic SNVs and short indels in the genes marked as cancer genes, matched by exact "
                 "position and allele, from the VCF.</li>"
                 "<li><b>Not assessed</b>: exon-level deletions/duplications and structural variants (a common mechanism for BRCA1, PMS2 and MSH2/EPCAM), "
                 "deep-intronic and regulatory variants, variants of uncertain significance (not interpreted as risk), mosaicism, and gene-level callable "
                 "coverage. PMS2, and other genes with pseudogenes or segmental duplications, are poorly resolved by short reads.</li>"
                 "<li><b>Absent sites</b> in a variant-only VCF are assumed reference. Confirm anything that matters in a clinical laboratory.</li></ul>")
    return parts


def _clinical_tests_html(res, geno, cv_hits) -> list:
    from . import clinical_tests
    rows = clinical_tests.summarize(res, geno, cv_hits)
    cls = {"targeted": "act-ok", "partial": "act-dose", "limited": "act-dose", "none": "act-bad"}
    body = "".join(f"<tr><td>{E(r['area'])}</td><td>{E(r['test'])}</td><td>{E(r['result'])}</td>"
                   f"<td><b class='{cls.get(r['coverage'], '')}'>{E(r['coverage'])}</b></td><td>{E(r['limits'])}</td></tr>" for r in rows)
    return ["<h2>Common clinical genetic tests: what this report covers</h2><p>A map from familiar test names to what was checked in this genome. "
            "<b>This is a research screen, not the test.</b> Coverage: <i>targeted</i> = the named variants are read directly; <i>partial</i> = ClinVar pathogenic/likely-pathogenic SNVs and short indels (2+ review stars) "
            "plus whatever optional modules were run, but not structural variants or variants of uncertain significance; <i>limited</i> = weak for technical reasons; <i>none</i> = cannot be derived from a variant file.</p>"
            "<table><tr><th>Area</th><th>Test</th><th>Result in this genome</th><th>Coverage</th><th>Limits</th></tr>" + body + "</table>"]


def _carrier_html(res: dict, cv_hits) -> list:
    """Carrier status for family planning (carrier.py). Common carrier alleles are not frequency-filtered."""
    import pandas as pd
    from . import carrier
    out = ["<h2>4b. Carrier status for family planning</h2>",
           "<p>Carrier status matters for a pregnancy when <b>both</b> genetic parents carry a variant in the same recessive gene, or when the "
           "person providing the X chromosome carries an X-linked one. Compare gene by gene; the other parent needs testing that also covers what "
           "this report cannot see. <b>This is not a laboratory carrier screen, and finding nothing for a gene does not rule out being a carrier.</b></p>"]
    if cv_hits is None:
        return out + ["<p>The ClinVar screen was skipped, so carrier status was not assessed.</p>"]
    found = carrier.findings(cv_hits, res.get("sex"))

    def table(rows, extra=()):
        df = pd.DataFrame(rows)
        df["inheritance"] = [r["moi"] + (" (also has a dominant condition: one copy may matter for your own health)" if r["also_dominant"] else "")
                             for r in rows]
        df["recessive / X-linked condition (gene level)"] = [_short(r["diseases"].replace("; ", "|"), 3) for r in rows]
        df["ClinVar condition for this variant"] = [_short(r["condition"], 3) for r in rows]
        df["ACMG carrier list"] = ["yes" if r["acmg_panel"] else "" for r in rows]
        return _table(df, ["gene", "inheritance", "ACMG carrier list", "rsid", "status", "ClinVar condition for this variant",
                           "recessive / X-linked condition (gene level)", "stars", *extra])

    out.append(f"<h3>Variants found in recessive or X-linked genes ({len(found['carrier'])})</h3>")
    out.append(table(found["carrier"]) if found["carrier"] else "<p>No known pathogenic or likely pathogenic SNV or short indel found.</p>")
    if found["carrier"]:
        out.append("<p class='mut'>The gene-level column lists what ClinGen curates for the gene; the ClinVar column says what this variant was "
                   "classified for. Check that they refer to the same condition.</p>")
    if found["risk_allele"]:
        out.append("<h4>Risk or low-penetrance alleles in these genes</h4><p class='mut'>Not simple carrier status: these raise the chance of a "
                   "condition by themselves, and may also matter in combination with another variant in the same gene.</p>" + table(found["risk_allele"]))
    low = res.get("carrier_low_review") or []
    low_found = carrier.findings(pd.DataFrame(low), res.get("sex")) if low else None
    if low_found and (low_found["carrier"] or low_found["risk_allele"]):
        out.append("<h4>Lower review status (fewer ClinVar review stars than the main screen)</h4><p class='mut'>Weaker evidence: leads to check, "
                   "not findings.</p>" + table(low_found["carrier"] + low_found["risk_allele"]))
    if found["uncurated"]:
        genes = ", ".join(sorted({h["gene"] for h in found["uncurated"]}))
        out.append(f"<p class='mut'>ClinVar hits in genes with no recessive curation here (see section 4; absence of a curation is not evidence "
                   f"of dominance): {E(genes)}.</p>")
    gaps = "".join(f"<tr><td>{E(g)}</td><td>{E(c)}</td><td>{E(s)}</td><td>{E(why)}</td></tr>" for g, (c, s, why) in carrier.GAPS.items() if g != "GBA")
    out.append("<h3>What this cannot tell you</h3><p>For these common carrier conditions the usual cause is partly or wholly invisible to this "
               "pipeline. If you are planning children, assume they are <b>not covered</b>.</p><table><tr><th>Gene</th><th>Condition</th>"
               f"<th>Status here</th><th>Why</th></tr>{gaps}</table>"
               "<p class='mut'>For every other gene too: exon or whole-gene deletions and duplications, rearrangements, most repeat changes, "
               "variants not yet in ClinVar, and new classifications after the release used are not assessed. Sites absent from a variant-only "
               "VCF are assumed to match the reference.</p>")
    rp = res.get("repeats")
    if rp:
        rows = [x for x in rp.get("loci", []) if x.get("interpreted") and (x["gene"] in ("FMR1", "FXN") or x["worst"] in ("intermediate range", "pathogenic range"))]
        body = "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['disease'])} ({E(x['inheritance'])})</td><td>{E(x['genotype'])}</td><td>{E(x['worst'])}</td>"
                       f"<td>{E(_repeat_note(x['gene'], x['worst'])) if x['worst'] != 'normal range' else ''}</td></tr>" for x in rows)
        out.append("<h3>Repeat lengths relevant to children</h3><p class='mut'>Fragile X and Friedreich ataxia, plus any other screened repeat "
                   "outside the normal range, including dominant ones whose length can change between generations.</p><table><tr><th>Gene</th>"
                   f"<th>Condition</th><th>Repeat counts</th><th>Category</th><th>Note</th></tr>{body}</table>")
    else:
        out.append("<p><b>Repeat expansions were not run</b> (fragile X and Friedreich ataxia repeats are not assessed).</p>")
    cov = carrier.panel_coverage(cv_hits)
    body = "".join(f"<tr><td>{E(r['gene'])}</td><td>{E(r['condition'])}</td><td>{E(r['result'])}</td></tr>" for r in cov)
    out.append(f"<details><summary>Gene-by-gene status for the {len(cov)} genes of the ACMG carrier-screening list (Gregg 2021)</summary>"
               f"<table><tr><th>Gene</th><th>Condition</th><th>What this report can say</th></tr>{body}</table></details>")
    out.append(f"<p class='mut'>Genes treated as recessive or X-linked: {len(carrier.recessive_genes()):,} (ClinGen Definitive, Strong or Moderate "
               "curations, plus the ACMG list). Recessive genes not yet curated by ClinGen appear only in section 4.</p>")
    return out


def _cnv_screen_html(nc: dict) -> str:
    from .cnv_screen import findings
    head = "<h2>7b. Actionable large deletions and duplications (read depth)</h2>"
    if nc.get("error"):
        return head + f"<p>{E(nc['error'])}</p>"
    found = findings(nc)
    flagged = [f for f in found if f["flag"]]
    unread = [f["name"] for f in found if f["call"] == "not assessable"]
    summary = ""
    if not flagged:
        summary = (f"<p><b>None of the {len(found) - len(unread)} assessed changes was detected.</b>"
                   + (f" Not assessable: {E('; '.join(unread))}." if unread else "") + "</p>")
    notes = []
    for f in flagged:
        text = f"<b>{E(f['name'])}: {E(f['call'])}</b> (associated with {E(f['condition'])})."
        if f["partial"]:
            text += " Only part of the region is affected: the information below applies to the whole recurrent change and may not apply."
        text += (f" <br><b>Documented follow-up</b> ({E(f['guidance'])}; {E(f['source'])}): {E(f['follow_up'])} "
                 "Summarised from excerpts: read the source.")
        if f["diseases"]:
            text += "<br>Other disease links: " + "; ".join(
                f"<b>{E(d['condition'])}</b> ({E(d['effect'])}; {E(d['evidence'])}; {E(d['source'])})" for d in f["diseases"]) + "."
        if f["psychiatric"]:
            text += "<br>" + E(f["psychiatric"])
        notes.append(f"<div class='warn'>{text}<br>This is a read-depth screen: confirm with chromosomal microarray or MLPA before acting.</div>")
    rows = "".join(f"<tr><td>{E(f['name'])}</td><td>{E({'del': 'deletion', 'dup': 'duplication'}[f['change']])}</td><td>{E(f['condition'])}</td>"
                   f"<td>{E(f['guidance'])}</td><td>{'&ndash;' if f['copies'] is None else format(f['copies'], '.2f')}</td><td>{E(f['call'])}</td></tr>"
                   for f in found)
    x = nc.get("x_control_copies")
    control = "" if x is None else (f" Built-in check: non-pseudoautosomal chrX reads {x:.2f} copies in this genome (expected about 1 with one X, "
                                    "about 2 with two).")
    return (head + "<p>Only changes with strong evidence, a substantial effect on a carrier and a documented follow-up are reported; the "
            "strength of that guidance differs and is shown. <b>This is not a genome-wide CNV search</b>: other deletions and duplications are "
            "not looked for.</p>" + summary + "".join(notes) +
            "<table><tr><th>Region</th><th>Change</th><th>Condition</th><th>Follow-up guidance</th><th>Median copies</th><th>Result</th></tr>"
            + rows + "</table>"
            f"<p class='mut'>Depth in 10 kb bins at MAPQ&ge;30 against an autosomal baseline of {nc['baseline_depth']:.0f}x." + control +
            " The method has detected known duplications in public genomes; it has not been tested on a known deletion. It can miss events under "
            "about 200 kb, mosaic events and exon-level changes in disease genes.</p>")


def write_report(out: Path, res: dict, results, dist_data, geno, cv_hits, panel_pcs, meta, user_pcs, reg) -> None:
    anc = res["ancestry"]
    parts = [f"<h1>Genome report</h1><p class='mut'>Input: {E(res['vcf'])} &middot; sex: {E(res['sex'])} "
             f"({E(res['sex_source'])}) &middot; generated by wgs-report</p>"]
    parts.insert(0, disclaimer.html())
    parts.append("<p class='mut'>Polygenic scores describe populations: even for height (the best case) a score explains about 40% of variance; for most traits much less. "
                 "Findings in the ClinVar section are research-grade screens, not clinical results.</p>")

    # --- ancestry
    fr = ", ".join(f"{k} {100 * v:.0f}%" for k, v in anc["superpop_fractions"].items())
    near = ", ".join(f"{k} ({v})" for k, v in list(anc["nearest_pops"].items())[:6])
    parts.append("<h2>1. Genetic ancestry</h2>")
    parts.append(f"<p>Nearest 1000 Genomes neighbours in PCA space: <b>{E(fr)}</b>; closest populations: {E(near)}. "
                 f"Distance from the {E(anc['superpop'])} centroid: {anc['distance_from_centroid_sd']:.1f} SD.</p>")
    if not res["ancestry_ok"]:
        parts.append("<div class='warn'>Your genome is not a clean match to a single reference group. The scores "
                     "below were trained on European-ancestry data and are less accurate, and possibly biased, "
                     "for other or admixed ancestries.</div>")
    elif anc["superpop"] != "EUR":
        parts.append("<div class='warn'>The scores below were trained mostly on European-ancestry data; their accuracy "
                     f"in {E(anc['superpop'])}-ancestry individuals is lower than the published values shown.</div>")
    parts.append(f"<img src='{_pca_fig(panel_pcs, meta, user_pcs)}' alt='PCA'>")
    if res.get("admixture"):
        g, e = res["admixture"].get("global", {}), res["admixture"].get("europe", {})
        parts.append("<h3>Supervised ADMIXTURE (1000 Genomes reference groups)</h3><p>Global (4 groups): " +
                     E(", ".join(f"{k} {100 * v:.1f}%" for k, v in g.items())) + ". Within Europe (5 populations): " +
                     E(", ".join(f"{k} {100 * v:.1f}%" for k, v in e.items())) + ". These are model-based proportions "
                     "against a handful of reference populations, not statements about where your ancestors lived. "
                     "Small percentages (a few %) are within the noise of the method, and neighbouring populations "
                     "(CEU vs GBR) are not separable.</p>")
    if res.get("mt_haplogroup"):
        parts.append(f"<p>Maternal (mtDNA) haplogroup: <b>{E(res['mt_haplogroup'])}</b> (haplogrep, PhyloTree 17.3).</p>")
    yh = res.get("y_haplogroup")
    if yh:
        parts.append(f"<p>Paternal (Y-chromosome) haplogroup: <b>{E(yh['haplogroup'])}</b> (Yleaf, defining marker {E(yh['marker'])}; {yh['valid_markers']:,} "
                     f"markers read, quality score {yh['qc']:.2f}). It follows the direct paternal line only and has no health meaning.</p>")
    elif res.get("mt_haplogroup"):
        parts.append("<p class='mut'>The Y-chromosome haplogroup has not been assigned.</p>")

    # --- summary table
    ok = [r for r in results if not r.get("skipped")]
    parts.append("<h2>2. Trait predictions from polygenic scores</h2>")
    parts.append(f"<p>Scores are standardized against {res['n_reference']} 1000 Genomes {E(res['reference_group'])} "
                 "individuals after removing the effect of the first 4 ancestry PCs. With ~500 reference people, percentiles in the "
                 "tails (below ~5 or above ~95) are imprecise. &ldquo;Typical for similar scores&rdquo; multiplies that z-score by the published validation "
                 "correlation: where people with a similar score fall on average, not a calibrated prediction for an individual. "
                 "<b>Precedence:</b> a measured value outranks a large-effect or Mendelian variant, which outranks the polygenic score. Only traits that pass the "
                 "heritability and prediction-accuracy rules are shown (rules and evidence: <code>traits.yaml</code>, <code>docs/METHODS.md</code>). "
                 "<b>How to read the score:</b> a <b>positive</b> score means a genetic tendency toward a <b>higher</b> value of that trait "
                 "(for diseases, a higher risk); a <b>negative</b> score means a tendency toward a <b>lower</b> value (lower risk). It describes "
                 "the polygenic score only, relative to people of your ancestry. It can disagree with a measured value or with a single "
                 "large-effect variant you carry, which the score does not capture (see the notes below each table).</p>")
    def _reading(r):
        name = r["label"].split(" (")[0]
        name = name if name[:2].isupper() else name[:1].lower() + name[1:]
        if abs(r["z"]) < 0.25:
            return f"about average {name}" if r["kind"] == "quantitative" else f"about average risk of {name}"
        d = "higher" if r["z"] > 0 else "lower"
        return f"{d} than average {name}" if r["kind"] == "quantitative" else f"{d} than average risk of {name}"

    def rows_for(lst):
        rows = ["<table><tr><th>Trait</th><th>Score (SD)</th><th>Percentile</th><th>Genetic tendency</th>"
                "<th>Typical for similar scores</th><th>Accuracy</th></tr>"]
        for r in lst:
            if r["kind"] == "quantitative":
                exp = (f"{r['abs_expected']:.1f} {E(r['abs_unit'])}" if "abs_expected" in r
                       else ("percentile only" if r["tier"] == "C" else f"{r['expected_sd']:+.2f} SD"))
                acc = f"R&sup2;&asymp;{r['r2']:.2f}"
            else:
                exp = "percentile only" if r["tier"] == "C" else f"odds &times;{r['relative_odds']:.2f}"
                acc = f"OR/SD {r['or_per_sd']}"
            gi0 = geno.set_index("rsid")
            eff = r.get("related_effect") or {}
            override = bool(r.get("related_not_in_score") and any(
                rs in gi0.index and (gi0.loc[rs, "n_effect"] or 0) > 0 and
                ((eff.get(rs, "up") == "up" and r["z"] < 0.25) or (eff.get(rs) == "down" and r["z"] > -0.25))
                for rs in (r.get("related_rsids") or [])))
            reading = E(_reading(r)) + (" <b>(you carry a large-effect variant the score does not include; see note)</b>" if override else "")
            carried_rel = [rs for rs in (r.get("related_rsids") or []) if rs in gi0.index and (gi0.loc[rs, "n_effect"] or 0) > 0]
            vr = variant_reading(r, carried_rel)
            if vr:
                reading = vr[0] + f" <span class='mut'>(genome-wide background: {E(_reading(r))})</span>"
                exp = vr[1]
            rows.append(f"<tr><td>{E(r['label'])}</td><td>{r['z']:+.2f}</td>"
                        f"<td>{_fmt_pct(r['percentile'])}</td><td>{reading}</td><td>{exp}</td><td>{acc}</td></tr>")
        rows.append("</table>")
        notes = []
        gi = geno.set_index("rsid")
        for r in lst:
            carried = [rs for rs in (r.get("related_rsids") or []) if rs in gi.index and (gi.loc[rs, "n_effect"] or 0) > 0]
            if carried:
                items = "; ".join(f"{E(rs)} {E(gi.loc[rs, 'genotype'])}: {E((r.get('related_text') or {}).get(rs) or gi.loc[rs, 'interpretation'])}" for rs in carried)
                notes.append(f"<p class='mut'><b>{E(r['label'])}:</b> the polygenic score reads &ldquo;{E(_reading(r))}&rdquo;, but you carry "
                             f"{items}. Single large-effect variants (and repeat-number effects the score cannot see) are not "
                             "captured by the score, so the measured value or your actual risk can differ."
                             + (" For this trait the carried variant(s) are <b>not among the score&rsquo;s variants</b>, so the score does not include their effect: "
                                "read them alongside it, not combined into one number." if r.get("related_not_in_score") else "") + "</p>")
        return "\n".join(rows + notes)

    def cards_for(lst):
        out_ = ["<div class='grid'>"]
        for r in lst:
            d = dist_data[r["id"]]
            out_.append(_trait_card(r, _dist_fig(d[0], d[1], d[2], d[3], d[4], r["label"])))
        out_.append("</div>")
        return "\n".join(out_)

    TIER_TITLES = {"A": "Tier A: strongest predictive performance (held-out R&sup2; &ge; 0.20, or strong binary performance)",
                   "B": "Tier B: moderate predictive performance (held-out R&sup2; 0.09&ndash;0.20)",
                   "C": "Tier C: weak prediction, percentile only (shown because it was requested)"}   # tiers rank published predictive strength only, not replication or evidence quality

    def by_tier(lst):
        out_ = []
        for t in ("A", "B", "C"):
            sub_ = [r for r in lst if r["tier"] == t]
            if sub_:
                out_.append(f"<h3>{TIER_TITLES[t]}</h3>")
                out_.append(rows_for(sub_))
                out_.append(cards_for(sub_))
        return "\n".join(out_)

    main = [r for r in ok if r.get("category") != "blood"]
    blood = [r for r in ok if r.get("category") == "blood"]
    parts.append(by_tier(main))
    for r in results:
        if r.get("skipped"):
            parts.append(f"<p class='mut'>{E(r['label'])}: not reported ({E(r['skipped'])}).</p>")
    if blood:
        parts.append("<h2>2b. Genetic tendencies for blood-test values</h2>"
                     "<p>These are <b>inherited tendencies</b>, not predictions of your blood test result. Each score "
                     "explains roughly 10&ndash;35% of the variation in the lab value in UK Biobank; what you actually measure "
                     "also depends on iron, B12/folate, inflammation, medication, hydration, thalassaemia traits and illness. "
                     "Read a percentile of 80 for MCV as &ldquo;genetically inclined toward larger red cells than average&rdquo;. "
                     "A real result outside the lab range needs a clinician, whatever this page says.</p>")
        parts.append(by_tier(blood))

    if res.get("measured"):
        parts.append("<h2>2c. Your measured results vs the genetic tendency</h2><p>Values from your own lab spreadsheet. "
                     "&ldquo;Approx. SD&rdquo; places the mean of your results on a rough adult-male scale (approximate "
                     "typical values, <i>not</i> your lab's reference range; analyser-dependent or skewed tests are compared "
                     "qualitatively). &ldquo;Genetic shift&rdquo; is the score's expected deviation in SD. Disagreement is "
                     "normal: diet, medication, iron, illness and measurement variation matter more than genetics for most "
                     "of these tests.</p><table><tr><th>Test</th><th>Your values</th><th>Approx. SD</th><th>Genetic percentile</th>"
                     "<th>Genetic shift (SD)</th><th>Reading</th></tr>")
        for m in res["measured"]:
            vals = ", ".join(f"{v:g}" for _, v in m["values"])
            if m.get("n_on_treatment"):
                vals += f" ({m['n_on_treatment']} on treatment, excluded from comparison)"
            z = "n/a" if m["approx_z_measured"] is None else f"{m['approx_z_measured']:+.1f}"
            parts.append(f"<tr><td>{E(m['label'])}</td><td>{E(vals)} {E(m['unit'])}</td><td>{z}</td>"
                         f"<td>{_fmt_pct(m['genetic_percentile'])}</td><td>{m['expected_sd']:+.2f}</td><td>{E(m['verdict'])}</td></tr>")
            if m.get("clinical_flag"):
                parts.append(f"<tr><td colspan='6'><b>{E(m['label'])}:</b> {E(m['clinical_flag'])}</td></tr>")
            if m.get("note"):
                parts.append(f"<tr><td colspan='6' class='mut'>{E(m['label'])}: {E(m['note'])}</td></tr>")
        parts.append("</table>")

    # --- top variants
    parts.append("<h3>Largest single-variant contributions</h3><p class='mut'>Variants whose genotype moves your score "
                 "furthest from the reference average (in SD of the score). Individual variants are tiny effects; "
                 "this is shown for curiosity, not interpretation.</p>")
    for r in ok:
        tv = r.get("top_variants") or []
        if not tv:
            continue
        parts.append(f"<details><summary>{E(r['label'])}</summary><table><tr><th>Position (GRCh38)</th><th>Ref&gt;Alt</th>"
                     "<th>Alt copies</th><th>Contribution (SD)</th></tr>" +
                     "".join(f"<tr><td>chr{E(str(v['chrom']))}:{v['pos']}</td><td>{E(v['ref'])}&gt;{E(v['alt'])}</td>"
                             f"<td>{v['genotype_alt_count']:.0f}</td><td>{v['contribution_sd']:+.3f}</td></tr>" for v in tv) +
                     "</table></details>")

    # --- SNP panel
    parts.append("<h2>3. Notable single SNPs</h2>"
                 "<p class='mut'>Evidence level: <b>A</b> clinically established and guideline-backed (but a screen of selected alleles only); "
                 "<b>B</b> robust association or well-established genotype&ndash;phenotype link with variable effect size or penetrance; "
                 "<b>C</b> small or inconsistent effect, interest only. Most single-SNP trait claims are small-effect, ancestry-dependent "
                 "associations from population studies: tendencies, not outcomes.</p>")
    for cat, title in (("trait", "Appearance and traits"), ("nutrition", "Nutrition and lifestyle"),
                       ("health", "Health-related markers"), ("pgx", "Pharmacogenomics (drug response)"),
                       ("sensitive", "Sensitive (APOE)")):
        g = geno[geno.category == cat]
        if len(g) == 0:
            continue
        parts.append(f"<h3>{title}</h3>")
        if cat == "sensitive" and res.get("apoe"):
            parts.append(f"<p>APOE genotype: <b>{E(res['apoe'])}</b></p>")
        parts.append(_table(g, ["rsid", "gene", "trait", "evidence", "genotype", "interpretation", "status"],
                            ["rsID", "Gene", "Trait", "Evidence", "Genotype", "Meaning", "Call"]))
    parts.append("<p class='mut'>Single-SNP tags cover only a few common alleles. HLA proxies are not a substitute for "
                 "HLA typing.</p>")
    parts.extend(_hla_html(res.get("hla")))
    if res.get("pgx"):
        parts.append("<h3>Pharmacogenomic summary (CPIC-style, simplified)</h3><table><tr><th>Gene</th><th>Result</th>"
                     "<th>Phenotype</th><th>Drugs affected</th><th>Type of test</th><th>Evidence</th></tr>" +
                     "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['result'])}</td><td>{E(x['phenotype'])}</td>"
                             f"<td>{E(x['drugs'])}</td><td>{E(x.get('test_type', ''))}</td><td>{E(x['evidence'])}</td></tr>" for x in res["pgx"]) +
                     "</table>" + _pgx_consequence_html(res["pgx"]) + "<p class='mut'>Limited screens look at a few common variants only: &ldquo;not detected&rdquo; is not the same as "
                     "&ldquo;normal&rdquo;, and untested variants leave residual risk. Research-grade results: do not change any medication on the strength of this page.</p>")

    parts.extend(_pharmcat_html(res.get("pharmcat")))
    ab = res.get("abo")
    if ab and ab.get("group") != "undetermined":
        rh = (res.get("copynumber") or {}).get("loci", {}).get("RHD + RHCE", {}).get("estimated_total_copies")
        rhtxt = "" if rh is None else (" Rh(D) from RHD copy number: estimated total RHD+RHCE copies " f"{rh:.1f}" +
                                       (" (consistent with Rh-positive)" if rh >= 3.5 else
                                        " (one RHD deleted, Rh-positive carrier of the deletion?)" if rh >= 2.5 else
                                        " (would fit Rh-negative; confirm by serology)") + ".")
        parts.append(f"<h3>Blood group (inferred)</h3><p>ABO genotype <b>{E(ab['genotype'])}</b>, i.e. blood group <b>{E(ab['group'])}</b>. "
                     f"{E(ab['disclaimer'])}{rhtxt}</p>")

    parts.extend(_qc_html(res.get("qc")))

    parts.extend(_repeats_html(res.get("repeats")))

    # --- ClinVar
    parts.append("<h2>4. Rare diseases and carrier status (ClinVar pathogenic / likely-pathogenic variants)</h2>")
    if cv_hits is None or len(cv_hits) == 0:
        parts.append("<p>No matching variants (or ClinVar screen skipped).</p>")
    else:
        h = cv_hits.sort_values(["acmg_sf_gene", "stars"], ascending=False).copy()
        h["condition"] = h["condition"].map(lambda x: _short(x, 3))
        h["type"] = h["significance"].map(lambda x: "risk allele or low-penetrance variant (not simple carrier status)"
                                         if ("risk" in str(x).lower() or "low_penetrance" in str(x).lower())
                                         else "pathogenic / likely pathogenic")
        h["significance"] = h["significance"].map(lambda x: _short(x, 2, "/"))
        a = h[h.acmg_sf_gene]
        parts.append("<h3>In ACMG secondary-findings genes</h3>")
        if len(a):
            parts.append(_table(a, ["gene", "rsid", "zygosity", "condition", "significance", "type", "stars", "consequence"]))
        else:
            parts.append("<p>None found by this research pipeline. This is <b>not</b> a negative clinical result: only ClinVar "
                         "pathogenic/likely-pathogenic variants with at least 2 review stars, matched by exact position and allele "
                         "(SNVs and short indels), are screened. Structural variants, copy-number changes, novel variants and variants "
                         "classified after the ClinVar release used are not assessed.</p>")
        b = h[~h.acmg_sf_gene]
        b = b[(b.af_1000g.isna()) | (b.af_1000g < 0.01)]
        b = b.assign(carrier_panel=b.get("carrier_panel_gene", False).map(lambda x: "yes" if x else ""))
        if "cancer_gene" not in b:
            b = b.assign(cancer_gene="")
        b = b.assign(cancer_gene=b["cancer_gene"].fillna("").map(lambda x: f"cancer gene ({x})" if x else ""))
        cols = ["gene", "carrier_panel", "cancer_gene", "rsid", "zygosity", "pattern", "condition", "significance", "type", "stars"]
        parts.append("<h3>Rare-disease findings in other genes (1000 Genomes AF &lt; 1%)</h3>")
        if "pattern" in b:
            hom = b[b.pattern == "homozygous"]
            cmp_ = b[b.pattern.str.startswith("possible compound")]
            sing = b[b.pattern == "single heterozygous"]
            for title, sub, note in (
                    ("Homozygous", hom, "Two copies of a pathogenic variant: if the condition is recessive this is the affected genotype, so review it first."),
                    ("Possible compound heterozygotes", cmp_, "Two different pathogenic variants in one gene. They may be on the same chromosome (no disease) or on opposite ones (affected if recessive); phase needs parents or long reads."),
                    ("Single heterozygous variants", sing, "Genes marked <i>carrier_panel</i> are on the ACMG tier-3 carrier-screening list (Gregg 2021; 113 genes), where one copy is carrier status for a recessive or X-linked condition. Genes marked <i>cancer_gene</i> have Definitive or Strong ClinGen evidence for cancer predisposition (the inheritance modes are listed; dominant means one copy raises risk, recessive means carrier status). Carrier status if the condition is recessive (relevant for family planning); possible disease risk if dominant.")):
                parts.append(f"<h4>{title} ({len(sub)})</h4><p class='mut'>{note}</p>" + (_table(sub, cols) if len(sub) else "<p>None.</p>"))
        else:
            parts.append(_table(b, cols))
        parts.append("<p class='mut'>How to read these: one copy of a variant in a <b>recessive</b>-disease gene usually means "
                     "carrier status (relevant for family planning, not for your own health). A <b>risk allele</b> (for example a "
                     "thrombophilia variant) is a different thing: it raises the chance of a condition without being a simple carrier "
                     "state, and the SNP section explains its relative risk. ClinVar classifications can change and are condition-specific; "
                     "inheritance is not assessed in this section (see 4b for carrier status, which also includes variants with 1000 Genomes AF of 1% or more, "
                     "hidden here). Confirm anything that matters in a clinical laboratory. Full table: "
                     "clinvar_pathogenic_hits.tsv.</p>")

    parts.extend(_cancer_html(res, cv_hits))
    parts.extend(_carrier_html(res, cv_hits))
    parts.extend(_clinical_tests_html(res, geno, cv_hits))

    # --- rare LoF, mtDNA
    if res.get("lof_scan_run"):
        parts.append("<h2>5. Rare predicted loss-of-function variants in medically relevant genes (not in ClinVar P/LP)</h2>"
                     "<p>Nonsense, frameshift, canonical-splice and start-loss variants in the ACMG actionable-gene list, the 113 ACMG carrier-screening genes, 119 ClinGen cancer-predisposition genes (Definitive or Strong evidence) and 21 genes where rare loss-of-function variants raise atrial-fibrillation risk (JAMA Cardiology 2025 supplement), "
                     "1000 Genomes AF &lt; 1%. A candidate here is <b>not</b> a clinically reportable variant. Before it means anything it needs "
                     "read-level inspection, allele-balance and repeat-context review, transcript review and orthogonal confirmation. Whole-genome data is clearly informative here, but these are "
                     "<b>candidates</b>: many ACMG genes cause disease by other mechanisms, some variants sit in "
                     "benign regions or are alignment artefacts at repeats.</p>")
        lc = res.get("lof_candidates") or []
        if lc:
            parts.append("<table><tr><th>Gene</th><th>Position (GRCh38)</th><th>Ref&gt;Alt</th><th>Genotype</th>"
                         "<th>Consequence</th><th>Gene list</th><th>Note</th></tr>" +
                         "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['chrom'])}:{x['pos']}</td><td>{E(x['ref'])}&gt;{E(x['alt'])}</td>"
                                 f"<td>{E(x['genotype'])}</td><td>{E(x['consequence'])}</td><td>{E(x.get('list') or '')}</td><td>{E(x.get('flag') or '')}</td></tr>"
                                 for x in lc) + "</table>")
        else:
            parts.append("<p>No candidates found.</p>")
    mt = res.get("mito")
    if mt:
        vs = mt["variants"]
        het = [v for v in vs if v["alt_fraction"] is not None and 0.1 < v["alt_fraction"] < 0.9]
        parts.append(f"<h2>6. Mitochondrial DNA</h2><p>Mean depth {mt['mean_depth']:.0f}x; {len(vs)} differences from the "
                     "rCRS reference (all listed variants are in the downloadable JSON). " +
                     (f"Haplogroup {E(res['mt_haplogroup'])} (see section 1). " if res.get("mt_haplogroup") else
                      "A haplogroup needs a dedicated tool such as haplogrep and is not assigned here. ") +
                     f"Possible heteroplasmy (alt fraction 10&ndash;90%): {len(het)} sites.</p>")

    cn = res.get("copynumber")
    if cn:
        from .copynumber import LOCI
        rows = "".join(f"<tr><td>{E(k)}</td><td>{v['estimated_total_copies']:.1f}</td><td>{E(LOCI.get(k, v)['expected'])}</td>"
                       f"<td>{E(LOCI.get(k, v).get('confidence', ''))}</td><td>{E(LOCI.get(k, v)['meaning'])}</td></tr>"
                       for k, v in cn["loci"].items())
        lp = cn["lpa"]
        parts.append("<h2>7. Repeat- and copy-number-dependent loci (read depth)</h2>"
                     "<p>Estimated from sequencing depth, with all reads counted so that near-identical copies are not "
                     "lost. Each value is a <b>diploid total</b> (both chromosomes) with an error of about &plusmn;1 copy at best. "
                     "More fundamentally, for near-identical genes (SMN1/SMN2, C4A/C4B, HBA1/HBA2) and repeat arrays (LPA KIV-2) the quantity "
                     "that matters &mdash; gene-specific dose, isotype, or the size of each chromosome&rsquo;s array &mdash; cannot be recovered "
                     "from total depth. Paralog-ratio assays, ddPCR, MLPA or long reads are needed; GC bias and subtelomeric mapping add "
                     "error.</p><table><tr><th>Locus</th><th>Estimated total copies</th>"
                     f"<th>Typical</th><th>Confidence</th><th>Interpretation and limits</th></tr>{rows}</table>"
                     f"<p><b>LPA / Lp(a):</b> read depth over the kringle IV repeat block is {lp['repeat_block_to_unique_ratio']:.2f}&times; "
                     "that of single-copy LPA exons. This ratio is not calibrated to a KIV-2 copy count. "
                     "Lp(a) is driven by the size of the <i>shorter</i> isoform, which a "
                     "total cannot reveal: a short chromosome with few repeats plus a long one can give the "
                     "same total as two intermediate ones. Treat this as context, not an explanation, for a measured Lp(a). "
                     "A value one copy away from typical is a hypothesis to confirm with a dedicated assay, not a finding.</p>")

    if cn and cn.get("cnv_screen"):
        parts.append(_cnv_screen_html(cn["cnv_screen"]))

    # --- appendix
    d = res["dosage_stats"]
    parts.append("<h2>Appendix. Data quality and method notes</h2><ul>"
                 f"<li>Variant records in VCF: {res['vcf_summary'].get('number of records', 'n/a'):,}; SNPs: "
                 f"{res['vcf_summary'].get('number of SNPs', 'n/a'):,}; ts/tv {res['vcf_summary'].get('ts/tv', 'n/a')}.</li>"
                 f"<li>Of {d['n_variants']:,} scored sites: {d['from_vcf_het']:,} heterozygous and {d['from_vcf_hom_alt']:,} "
                 f"homozygous-alt calls; {d['assumed_hom_ref']:,} absent from the VCF and treated as homozygous reference "
                 f"(the VCF lists variants only); {d['missing_filtered_or_nocall']:,} filtered/no-call and mean-imputed.</li>"
                 f"<li>Callability check (a variant-only VCF cannot prove absent sites are reference): hom-ref {100 * res['callability']['obs_homref']:.1f}% observed vs "
                 f"{100 * res['callability']['exp_homref']:.1f}% expected from 1000G allele frequencies; het {100 * res['callability']['obs_het']:.1f}% vs "
                 f"{100 * res['callability']['exp_het']:.1f}%; hom-alt {100 * res['callability']['obs_homalt']:.1f}% vs {100 * res['callability']['exp_homalt']:.1f}%. "
                 "This global check cannot show that any single site or region is callable: without a gVCF, reference-confidence blocks or a "
                 "callable-regions mask, &ldquo;absent from the VCF&rdquo; may mean homozygous reference, low depth, a masked region or a no-call. "
                 "Difficult regions (segmental duplications, HLA, centromeres) can therefore bias a polygenic score slightly; the curated SNP "
                 "panel is checked against the sequencing reads instead.</li>"
                 "<li>&ldquo;Typical for similar scores&rdquo; uses E[trait | score] = &radic;R&sup2; &times; z. No individual interval is shown because "
                 "an interval under a simplified bivariate-normal model would not be a calibrated personal prediction interval. "
                 "It would omit: error in the published R&sup2;; transferability to your ancestry; age and sex effects; differences in "
                 "phenotype and assay definitions; environment; calibration in an individual; within-family versus population-level accuracy; "
                 "selection in the UK Biobank; and uncertainty from the ancestry-PC adjustment.</li>"
                 "<li>Odds multipliers for disease scores are OR<sup>z</sup> relative to an average-score person. They are not absolute risks "
                 "and must not be read as a percent change in probability; an absolute risk needs the baseline incidence for your age, sex, "
                 "ancestry and country plus calibration.</li>"
                 "<li>Percentiles come from a Gaussian model fitted to about 500 reference people; the plausible range shown beside each one "
                 "reflects that sampling error only. Tail percentiles are especially imprecise.</li>"
                 "<li>R&sup2; values come from UK Biobank or similar cohorts; within-family (sibling) accuracy is lower, and accuracy "
                 "drops further outside European ancestry.</li></ul>")

    pv = res.get("provenance")
    if pv:
        tv = ", ".join(f"{k} {E(str(v))}" for k, v in (pv.get("tools") or {}).items() if v)
        parts.append(f"<p class='small'>Provenance: wgs-report {E(str(pv.get('wgsreport')))}, {E(str(pv.get('date')))}, Python {E(str(pv.get('python')))}, "
                     f"{tv}; {len(pv.get('pgs_score_sha256_16') or {})} PGS files (checksums in results.json)."
                     + (" <b>Reference-build check was skipped.</b>" if pv.get("build_check_skipped") else "") + "</p>")
    body = "\n".join(parts)
    (out / "report.html").write_text(
        f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>Genome report</title><style>{CSS}</style></head><body>{body}</body></html>")
