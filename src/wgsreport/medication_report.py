"""Printable drug-response summary (pharmacogenomics) for the person whose genome it is.

Self-contained HTML with print CSS. It carries no sample ID or name. The header says how the results were produced so the reader can judge them.
It is research output and does not suggest sharing it with anyone."""
from __future__ import annotations

import html
from pathlib import Path

from . import disclaimer
from .report import ACTION_CLASS

E = html.escape

CSS = """
body{font:14px/1.45 system-ui,Arial,sans-serif;color:#111;background:#fff;max-width:900px;margin:0 auto;padding:20px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:22px 0 6px;border-bottom:1px solid #999;padding-bottom:2px}
table{border-collapse:collapse;width:100%;margin:6px 0;font-size:13px}th,td{border:1px solid #bbb;padding:4px 6px;text-align:left;vertical-align:top}
th{background:#eee}.warn{border:1px solid #c90;background:#fff4d6;padding:8px 12px;margin:10px 0}.box{border:1px solid #888;padding:8px 12px;margin:10px 0;background:#fafafa}.mut{color:#555;font-size:12px}
.act-ok{color:#1a7f37}.act-dose{color:#9a6700}.act-bad{color:#cf222e}.form{display:flex;gap:24px;margin:10px 0}.form span{flex:1;border-bottom:1px solid #000;padding-top:18px}
@media print{body{padding:0;font-size:12px}h2{break-after:avoid}tr{break-inside:avoid}.noprint{display:none}}
"""

# Per gene: guideline source, what the screen does NOT cover, and the confirmatory clinical test.
GENE_INFO = {
    "CYP2D6": ("CPIC: codeine/tramadol (2021), tricyclics (2016), SSRIs (2023), tamoxifen (2018)",
               "structural variants, gene deletions/duplications and rare alleles depend on the star-allele caller", "clinical CYP2D6 genotyping with copy-number analysis"),
    "CYP2C19": ("CPIC: clopidogrel (2022), PPIs (2020), citalopram/escitalopram/sertraline (2023), voriconazole (2016)",
                "only *2, *3 and *17 are tested; rarer no-function alleles (*4-*8) are not", "clinical CYP2C19 genotype panel"),
    "CYP2C9": ("CPIC: warfarin (2017), phenytoin (2020), NSAIDs (2020)", "only *2 and *3 are tested; *5, *6, *8, *11 are not", "clinical CYP2C9 genotype panel"),
    "VKORC1": ("CPIC: warfarin (2017)", "one promoter variant; CYP4F2 and clinical factors also affect dose", "use a pharmacogenetic dosing algorithm with the clinical genotype"),
    "SLCO1B1": ("CPIC: statins (2022)", "only *5 (c.521T>C) is tested", "clinical SLCO1B1 genotype"),
    "DPYD": ("CPIC: fluoropyrimidines (2017, updated 2024)", "three variants only (*2A, D949V, *13); other reduced-function variants (e.g. HapB3) are not tested, so a normal result does not exclude DPD deficiency",
             "full DPYD genotype or DPD enzyme testing before treatment"),
    "TPMT": ("CPIC: thiopurines (2018)", "only *2, *3B, *3C are tested", "TPMT enzyme activity or full genotype"),
    "NUDT15": ("CPIC: thiopurines (2018)", "only *3 (R139C) is tested", "clinical NUDT15 genotype"),
    "G6PD": ("CPIC: rasburicase and other oxidant drugs (2022)", "one variant (A-); Mediterranean and many other deficiency variants are not tested", "quantitative G6PD enzyme assay"),
    "HLA-B*57:01": ("CPIC: abacavir (2020)", "a proxy SNP, not HLA typing", "clinical HLA-B*57:01 test"),
    "CYP3A5": ("CPIC: tacrolimus (2015)", "only *3 is tested; *6 and *7 (mostly African ancestry) are not", "clinical CYP3A5 genotype"),
    "CYP2B6": ("CPIC: efavirenz (2019)", "only *6 is tested; *18 and others are not", "clinical CYP2B6 genotype"),
    "UGT1A1": ("CPIC: atazanavir (2015); DPWG for irinotecan", "only *6 is tested; the common *28 promoter repeat is not tested", "clinical UGT1A1 genotype including *28"),
}


def _act(a: str) -> str:
    return f"<b class='{ACTION_CLASS.get(a, '')}'>{E(a)}</b>"


def _key(gene: str) -> str:
    for k in GENE_INFO:
        if gene.startswith(k.split(" ")[0].split("(")[0]) or gene.split(" ")[0] == k:
            return k
    return gene


def build(res: dict) -> str:
    pgx = res.get("pgx") or []
    pv = res.get("provenance") or {}
    pc = res.get("pharmcat")
    tools = ", ".join(f"{k} {v}" for k, v in (pv.get("tools") or {}).items() if v)
    flagged, routine = [], []
    for x in pgx:
        for c in x.get("consequences") or []:
            (routine if c["action"] == "use normally" else flagged).append((x, c))
    o = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
         f"<title>Drug response summary</title><style>{CSS}</style></head><body>",
         "<h1>Drug response summary</h1>", disclaimer.html(),
         "<p class='mut'>Research-grade result from whole-genome sequencing. Please read &ldquo;How this report was produced&rdquo; below.</p>",
         "<h2>How this report was produced</h2><div class='box'><ul>"
         f"<li><b>Source data:</b> a whole-genome sequencing variant file (VCF, GRCh38) from a consumer sequencing service, analysed locally with the open-source tool wgsreport {E(str(pv.get('wgsreport', '')))} "
         f"({E(str(pv.get('date', '')))}). Tools: {E(tools) or 'bcftools, samtools, plink2'}.</li>"
         "<li><b>Method:</b> specific variants that define common drug-metabolism alleles were read from the file and converted to a "
         "metabolizer phenotype and drug advice, simplified from <b>CPIC</b> guidelines (cpicpgx.org)"
         + (f"; star-allele calls and recommendations from <b>PharmCAT {E(str(pc.get('version') or ''))}</b> are included below" if pc else "") + ".</li>"
         "<li><b>Not a clinical test:</b> this was not run in an accredited (CLIA/ISO 15189) laboratory and is not a medical record. "
         "Only the listed variants were examined; &ldquo;no variant detected&rdquo; does not exclude other variants, copy-number or structural changes. "
         "Genes marked as screens (not full star-allele typing) are limited.</li>"
         "<li><b>Use:</b> for your own information only. Do not start, stop or change any medicine on the strength of this page.</li></ul></div>"]
    n_flag = len({x["gene"] for x, _ in flagged})
    o.append(f"<h2>Summary</h2><div class='box'><b>{len(pgx)}</b> drug-metabolism genes were screened. <b>{n_flag}</b> have a result that may change "
             "drug response: " + (", ".join(sorted({x["gene"] for x, _ in flagged})) or "none") + ". "
             "The tables are: (1) results that may change drug response, (2) genes with no change, (3) an alphabetical drug lookup, "
             "(4) gene-by-gene detail with what was <i>not</i> tested and which clinical test would confirm it"
             + (", (5) PharmCAT calls" if pc else "") + ".</div>")
    o.append("<h2>1. Results that may change drug response</h2>")
    if flagged:
        o.append("<table><tr><th>Gene</th><th>My result</th><th>Drug</th><th>Consequence</th><th>Usual advice</th></tr>" +
                 "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(x['result'])}<br><span class='mut'>{E(x['phenotype'])}</span></td><td>{E(c['drug'])}</td>"
                         f"<td>{E(c['effect'])}</td><td>{_act(c['action'])}</td></tr>" for x, c in flagged) + "</table>")
    else:
        o.append("<p>No tested variant changes the usual advice for the drugs covered.</p>")
    o.append("<h2>2. Genes where no change to usual drug response was found</h2>")
    if routine:
        o.append("<table><tr><th>Gene</th><th>Drugs</th><th>Note</th></tr>" +
                 "".join(f"<tr><td>{E(x['gene'])}</td><td>{E(c['drug'])}</td><td>{E(c['effect'])}</td></tr>" for x, c in routine) + "</table>")
    else:
        o.append("<p>None.</p>")
    # alphabetical drug lookup: one row per drug name in each consequence (comma-separated lists are split)
    look = {}
    for x in pgx:
        for c in x.get("consequences") or []:
            for d in [d.strip() for d in c["drug"].split(",")]:
                look.setdefault(d.lower(), []).append((x["gene"], c))
    o.append("<h2>3. Drug lookup (alphabetical)</h2><table><tr><th>Drug</th><th>Gene</th><th>Usual advice</th><th>Consequence</th></tr>" +
             "".join(f"<tr><td>{E(d)}</td><td>{E(g)}</td><td>{_act(c['action'])}</td><td>{E(c['effect'])}</td></tr>"
                     for d in sorted(look) for g, c in look[d]) + "</table>"
             "<p class='mut'>Drugs not listed here were not assessed; absence from the list does not mean no gene affects the drug.</p>")
    o.append("<h2>4. Gene-by-gene detail</h2><table><tr><th>Gene</th><th>Variants read (copies of the effect allele)</th><th>Result and phenotype</th>"
             "<th>Guideline</th><th>Not covered by this screen</th><th>Confirmatory test</th></tr>" +
             "".join("<tr><td>{g}</td><td>{v}</td><td>{r}<br><span class='mut'>{p}</span></td><td>{gl}</td><td>{nc}</td><td>{ct}</td></tr>".format(
                 g=E(x["gene"]), v=E(", ".join(f"{v['rsid']}: {v['copies']}" for v in x.get("variants") or []) or "star-allele caller"),
                 r=E(x["result"]), p=E(x["phenotype"]),
                 gl=E(GENE_INFO.get(_key(x["gene"]), ("", "", ""))[0]), nc=E(GENE_INFO.get(_key(x["gene"]), ("", "", ""))[1]),
                 ct=E(GENE_INFO.get(_key(x["gene"]), ("", "", ""))[2])) for x in pgx) + "</table>")
    if pc:
        o.append("<h2>5. PharmCAT star-allele calls</h2><table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th></tr>" +
                 "".join(f"<tr><td>{E(g['gene'])}</td><td>{E(g['diplotype'])}</td><td>{E(g['phenotype'])}</td></tr>"
                         for g in pc.get("genes", []) if not g["diplotype"].startswith("Unknown")) + "</table>")
        if pc.get("drugs"):
            o.append("<table><tr><th>Drug</th><th>Gene</th><th>CPIC recommendation</th><th>Strength</th></tr>" +
                     "".join(f"<tr><td>{E(d['drug'])}</td><td>{E(', '.join(d['genes']))}</td><td>{E(d['recommendation'])}</td><td>{E(d['strength'])}</td></tr>"
                             for d in pc["drugs"]) + "</table>")
    fl = (res.get("hla") or {}).get("flagged") or []
    if fl:
        o.append("<h2>HLA alleles with drug associations (typed from reads)</h2><table><tr><th>Allele</th><th>Copies</th><th>Meaning</th></tr>" +
                 "".join(f"<tr><td>{E(f['allele'])}</td><td>{f['copies']}</td><td>{E(f['meaning'])}</td></tr>" for f in fl) + "</table>")
    o.append("<h2>Notes</h2><ul><li>Allele frequencies differ by ancestry; variants not tested here (for example *28 in UGT1A1, "
             "HapB3 in DPYD, CYP2D6 duplications) may matter more in some populations.</li>"
             "<li>These are single-gene effects only. They do not account for interacting drugs (for example CYP2D6 or CYP2C19 inhibitors), "
             "kidney or liver function, age or other clinical factors.</li>"
             "<li>Where the row says &ldquo;confirm first&rdquo;, a clinical test is needed before relying on the result. "
             "&ldquo;Adjust dose&rdquo; and &ldquo;avoid&rdquo; summarise general guideline advice; they are not instructions.</li>"
             "<li>Result classes are the CPIC&nbsp;terms: &ldquo;use normally&rdquo; means no tested variant changes standard dosing.</li></ul>")
    o.append("<p class='mut'>Generated locally by wgsreport. No name or identifier is stored in this file. "
             "Guidelines: Clinical Pharmacogenetics Implementation Consortium (CPIC). Advice shown is general information.</p></body></html>")
    return "".join(o)


def write(out: Path, res: dict) -> Path:
    p = Path(out) / "medication_report.html"
    p.write_text(build(res))
    return p
