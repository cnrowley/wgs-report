"""Pharmacogenomic summary: genotype -> metabolizer phenotype -> drugs affected.

Simplified from CPIC/DPWG tables for the common alleles on the SNP panel (CYP2D6 comes from an external star-allele
caller such as Cyrius). Results assume no other star alleles are present; this is an educational summary, not a
clinical pharmacogenomic report.
"""
from __future__ import annotations

import pandas as pd

# CYP2D6 activity values (CPIC); only common alleles are mapped, anything else is reported as unmapped.
D6_ACTIVITY = {"*1": 1.0, "*2": 1.0, "*33": 1.0, "*35": 1.0, "*9": 0.5, "*17": 0.5, "*29": 0.5, "*41": 0.5,
               "*10": 0.25, "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0, "*7": 0.0, "*8": 0.0, "*11": 0.0,
               "*15": 0.0, "*36": 0.0, "*68": 0.0}
D6_DRUGS = ("codeine, tramadol (prodrugs; weaker effect if reduced), tamoxifen, metoprolol, atomoxetine, "
            "many antidepressants/antipsychotics (e.g. amitriptyline, paroxetine, aripiprazole)")

# Plain-language consequence of each result, per drug: (drug, what it means for you, action).
# Actions: USE = use normally; DOSE = dose adjustment usually advised; ALT = avoid, or choose an alternative; LESS = may work less well;
# CONFIRM = confirm with a definitive test first. Summarised from CPIC guidelines; general information only.
USE, DOSE, ALT, LESS, CONFIRM = ("use normally", "adjust dose", "avoid / choose alternative", "may be less effective", "confirm first")


def _c(drug, effect, action):
    return {"drug": drug, "effect": effect, "action": action}


def consequences_cyp2c19(loss: int, n17: int) -> list[dict]:
    if loss >= 2:
        return [_c("clopidogrel", "you activate it poorly, so it protects less against clots", ALT),
                _c("omeprazole, lansoprazole, pantoprazole", "you clear them slowly, so exposure is higher (more effect and more side effects)", DOSE),
                _c("citalopram, escitalopram", "you clear them slowly, so blood levels are higher", DOSE),
                _c("sertraline", "slower clearance; standard start, lower maintenance dose may suffice", DOSE),
                _c("voriconazole", "you clear it slowly, so levels are higher", ALT)]
    if loss == 1:
        return [_c("clopidogrel", "you activate it less well, so it protects less against clots", ALT),
                _c("omeprazole, lansoprazole, pantoprazole", "somewhat higher exposure; usual start, review dose on long-term use", USE),
                _c("citalopram, escitalopram, sertraline", "somewhat higher levels; usual start", USE),
                _c("voriconazole", "usual dose, with drug-level monitoring", USE)]
    if n17:
        return [_c("clopidogrel", "you activate it normally or slightly more", USE),
                _c("omeprazole, lansoprazole, pantoprazole", "you clear them faster, so they may work less well at the usual dose", LESS),
                _c("citalopram, escitalopram", "faster clearance; may be less effective, consider another drug if there is no response", LESS),
                _c("voriconazole", "faster clearance, so levels may be too low", LESS if n17 == 1 else ALT)]
    return [_c("clopidogrel, PPIs, citalopram/escitalopram, sertraline, voriconazole", "no tested change in metabolism", USE)]


def consequences_cyp2c9(score: float) -> list[dict]:
    if score >= 2:
        return [_c("warfarin, phenytoin, celecoxib, ibuprofen, flurbiprofen", "no tested reduced-function allele: usual metabolism", USE)]
    sev = "markedly" if score < 1 else "moderately"
    return [_c("warfarin", f"you clear it {sev} more slowly, so you need a lower dose (combine with VKORC1)", DOSE),
            _c("phenytoin", f"you clear it {sev} more slowly, so levels run higher and toxicity is more likely", DOSE),
            _c("celecoxib, flurbiprofen, ibuprofen, meloxicam, piroxicam", f"you clear them {sev} more slowly; start at a reduced dose", DOSE)]


def consequences_vkorc1(v: int) -> list[dict]:
    if v == 0:
        return [_c("warfarin", "usual sensitivity: usual dose requirement", USE)]
    return [_c("warfarin", "more sensitive to warfarin, so a lower dose is needed" + (" (strongly)" if v == 2 else ""), DOSE)]


def consequences_slco1b1(v: int) -> list[dict]:
    if v == 0:
        return [_c("simvastatin and other statins", "usual statin uptake into the liver", USE)]
    return [_c("simvastatin", "less is taken up by the liver, so blood levels and muscle-injury risk are higher" + (" (substantially)" if v == 2 else ""),
               DOSE if v == 1 else ALT),
            _c("atorvastatin, rosuvastatin, pravastatin", "smaller effect than with simvastatin; a lower dose or pravastatin/fluvastatin is often preferred", DOSE)]


def consequences_dpyd(n: int) -> list[dict]:
    if n == 0:
        return [_c("fluorouracil, capecitabine", "no tested reduced-function variant: usual dosing (untested variants remain possible)", USE)]
    if n == 1:
        return [_c("fluorouracil, capecitabine", "you break the drug down slowly, so severe or fatal toxicity is more likely at a full dose; "
                   "CPIC advises starting at about half the dose and adjusting", DOSE)]
    return [_c("fluorouracil, capecitabine", "very low breakdown capacity: serious toxicity risk even at reduced doses; avoid, or use only under specialist advice", ALT)]


def consequences_thiopurine(n: int, gene: str) -> list[dict]:
    drugs = "azathioprine, mercaptopurine, thioguanine"
    if n == 0:
        return [_c(drugs, f"no tested reduced-function {gene} allele: usual dosing", USE)]
    if n == 1:
        return [_c(drugs, "you inactivate it less well, so the bone marrow is more sensitive; CPIC advises a reduced starting dose (about 30-80%)", DOSE)]
    return [_c(drugs, "very little inactivation: severe bone-marrow suppression at usual doses; use a drastically reduced dose or an alternative", ALT)]


def consequences_g6pd(v: int) -> list[dict]:
    if v == 0:
        return [_c("rasburicase, primaquine, dapsone, nitrofurantoin", "A- variant not detected (other G6PD variants are not tested)", USE)]
    return [_c("rasburicase", "risk of red-cell destruction (haemolysis); contraindicated if G6PD-deficient", ALT),
            _c("primaquine, tafenoquine, dapsone, nitrofurantoin", "oxidant drugs that can cause haemolysis; avoid or use only with G6PD testing and monitoring", ALT)]


def consequences_hlab5701(v: int) -> list[dict]:
    if v == 0:
        return [_c("abacavir, flucloxacillin", "tag not detected; B*57:01 unlikely (HLA typing is definitive)", USE)]
    return [_c("abacavir", "if you carry HLA-B*57:01 you risk a serious hypersensitivity reaction; do not start until HLA typing is negative", CONFIRM),
            _c("flucloxacillin", "if you carry HLA-B*57:01, liver injury is more likely; confirm with HLA typing", CONFIRM)]


def consequences_cyp2d6(activity: float) -> list[dict]:
    if activity == 0:
        return [_c("codeine, tramadol", "you cannot convert them to the active form, so they give little or no pain relief", ALT),
                _c("amitriptyline, nortriptyline, paroxetine, fluvoxamine", "you clear them slowly, so levels are higher; lower dose or a different drug", DOSE),
                _c("tamoxifen", "less conversion to the active form, so possibly less effective; consider an alternative", LESS),
                _c("metoprolol, atomoxetine, aripiprazole", "higher levels; lower dose or slower titration", DOSE)]
    if activity < 1.25:
        return [_c("codeine, tramadol", "reduced conversion to the active form; may give less pain relief, switch if ineffective", LESS),
                _c("amitriptyline, nortriptyline, paroxetine", "somewhat higher levels; consider a lower starting dose", DOSE),
                _c("tamoxifen", "reduced conversion to the active form; CPIC suggests considering an alternative or avoiding CYP2D6 inhibitors", LESS)]
    if activity > 2.25:
        return [_c("codeine, tramadol", "you convert them too fast, so the risk of opioid toxicity is higher; avoid", ALT),
                _c("amitriptyline, nortriptyline, paroxetine", "lower levels; may be less effective", LESS)]
    return [_c("CYP2D6 drugs (codeine, tramadol, tricyclics, tamoxifen)", "usual metabolism by activity score", USE)]


def consequences_cyp3a5(n3: int) -> list[dict]:
    if n3 == 2:
        return [_c("tacrolimus", "you do not express CYP3A5 (*3/*3, the most common result in Europeans): usual starting dose", USE)]
    return [_c("tacrolimus", "you express CYP3A5, so you clear tacrolimus faster and the usual dose may give low levels; "
               "CPIC advises a higher starting dose (about 1.5-2x) with drug-level monitoring", DOSE)]


def consequences_cyp2b6(n6: int) -> list[dict]:
    if n6 == 0:
        return [_c("efavirenz, methadone", "no *6 allele: usual metabolism (other CYP2B6 alleles are not tested)", USE)]
    if n6 == 1:
        return [_c("efavirenz", "somewhat higher levels, so more central-nervous-system side effects; usual dose with monitoring", USE),
                _c("methadone", "somewhat slower clearance; titrate carefully", DOSE)]
    return [_c("efavirenz", "much higher levels and side effects; CPIC suggests a lower dose (400 mg/day) or an alternative", DOSE),
            _c("methadone", "slower clearance; titrate carefully", DOSE)]


def consequences_ugt1a1(n6: int) -> list[dict]:
    if n6 == 0:
        return [_c("irinotecan, atazanavir", "no *6 allele (the common *28 repeat variant is not tested here)", USE)]
    return [_c("irinotecan", "reduced inactivation of the active metabolite, so neutropenia and diarrhoea are more likely; consider a reduced dose", DOSE),
            _c("atazanavir", "higher bilirubin (jaundice) is more likely; usually harmless but may lead to stopping", USE)]


def parse_cyp2d6(genotype: str) -> dict:
    """Star-allele diplotype string ('hap1/hap2', '+' for tandems, 'xN' for copies) -> activity score and CPIC-style phenotype."""
    total, unmapped = 0.0, []
    for hap in genotype.split("/"):
        h = hap.strip()
        mult = 1
        if "x" in h:
            h, n = h.split("x")
            mult = int(n)
        parts = [p for p in h.split("+") if p]
        if len(parts) > 1:                        # tandem/hybrid arrangements: *68 (CYP2D7 hybrid) carries no function
            val = 0.0 if "*68" in parts else sum(D6_ACTIVITY.get(p, 0.0) for p in parts)
        else:
            val = D6_ACTIVITY.get(parts[0]) if parts else None
        if val is None:
            unmapped.append(h)
            val = 0.0
        total += val * mult
    if total == 0:
        ph = "poor metabolizer"
    elif total < 1.25:
        ph = "intermediate metabolizer"
    elif total <= 2.25:
        ph = "normal metabolizer"
    else:
        ph = "ultrarapid metabolizer"
    return {"activity_score": total, "phenotype": ph, "unmapped": unmapped}


def _cp(n: int) -> str:
    return "copy" if n == 1 else "copies"


def _n(g: pd.DataFrame, rsid: str):
    r = g[g.rsid == rsid]
    if len(r) == 0 or pd.isna(r.iloc[0]["n_effect"]):
        return None
    return int(r.iloc[0]["n_effect"])


def _src(g: pd.DataFrame, rsids: list[str]) -> str:
    s = {str(g[g.rsid == r].iloc[0]["status"]).split(" (")[0] for r in rsids if (g.rsid == r).any()}
    if s and all("confirmed from reads" in x for x in s):
        return "confirmed from reads"
    return ", ".join(sorted(s)) or "n/a"


def summarize(geno: pd.DataFrame, cyp2d6: str | None = None) -> list[dict]:
    rows = []

    NEG = "no tested reduced-function allele detected; residual risk from untested variants remains"

    def add(gene, result, phenotype, drugs, rsids=(), test_type="limited screen", conseq=()):
        rows.append({"gene": gene, "result": result, "phenotype": phenotype, "drugs": drugs, "test_type": test_type, "consequences": list(conseq),
                     "variants": [{"rsid": r, "copies": _n(geno, r)} for r in rsids],
                     "evidence": _src(geno, list(rsids)) if rsids else "star-allele caller"})

    if cyp2d6:
        d = parse_cyp2d6(cyp2d6)
        extra = f" (unmapped alleles: {', '.join(d['unmapped'])})" if d["unmapped"] else ""
        add("CYP2D6", cyp2d6, f"{d['phenotype']} (activity score {d['activity_score']:g}; classified as intermediate under the CPIC 2019 "
            f"consensus){extra}. Hybrid and other structural alleles should be confirmed by a clinical CYP2D6 assay before use", D6_DRUGS,
            test_type="structural call from a star-allele caller (needs specialised validation)",
            conseq=consequences_cyp2d6(d["activity_score"]))

    # CYP2C19: *2 rs4244285, *3 rs4986893 (no function); *17 rs12248560 (increased)
    n2, n3, n17 = _n(geno, "rs4244285"), _n(geno, "rs4986893"), _n(geno, "rs12248560")
    if None not in (n2, n3, n17):
        loss = n2 + n3
        if loss == 2:
            ph = "poor metabolizer"
        elif loss == 1:
            ph = "intermediate metabolizer" + (" (carries *17 on the other allele)" if n17 else "")
        elif n17 == 2:
            ph = "ultrarapid metabolizer"
        elif n17 == 1:
            ph = "rapid metabolizer"
        else:
            ph = "no tested reduced- or increased-function allele detected (normal metabolizer only if no untested alleles)"
        res = f"*2 x{n2}, *3 x{n3}, *17 x{n17}"
        add("CYP2C19", res, ph, "clopidogrel (needs activation; poor/intermediate respond less), PPIs (omeprazole etc.), "
            "citalopram/escitalopram, sertraline, voriconazole", ["rs4244285", "rs4986893", "rs12248560"],
            test_type="limited screen (*2, *3, *17 only)", conseq=consequences_cyp2c19(loss, n17))

    a, b = _n(geno, "rs1799853"), _n(geno, "rs1057910")
    if None not in (a, b):
        score = 2 - 0.5 * a - 1.0 * b
        ph = ("no tested reduced-function allele detected" if score >= 2 else
              "intermediate metabolizer" if score >= 1 else "poor metabolizer")
        add("CYP2C9", f"*2 x{a}, *3 x{b}", f"{ph} (activity score {score:g} from the alleles tested)",
            "warfarin, phenytoin, NSAIDs such as celecoxib/flurbiprofen", ["rs1799853", "rs1057910"],
            test_type="limited screen (*2, *3 only)", conseq=consequences_cyp2c9(score))
    v = _n(geno, "rs9923231")
    if v is not None:
        add("VKORC1 (-1639)", f"{v} {_cp(v)} of the lower-dose allele",
            ["usual warfarin sensitivity", "increased warfarin sensitivity (lower dose)",
             "high sensitivity (much lower dose)"][v], "warfarin", ["rs9923231"],
            test_type="single-SNP genotype (CYP4F2 and other factors also affect dose)", conseq=consequences_vkorc1(v))
    v = _n(geno, "rs4149056")
    if v is not None:
        add("SLCO1B1 *5", f"{v} {_cp(v)} of c.521C", ["no c.521T>C detected (this variant is the main known genetic risk factor for simvastatin myopathy)",
                                                     "decreased function (higher simvastatin myopathy risk)",
                                                     "poor function (avoid high-dose simvastatin)"][v],
            "simvastatin, atorvastatin and other statins", ["rs4149056"], test_type="single variant (*5 only)", conseq=consequences_slco1b1(v))
    v = _n(geno, "rs776746")
    if v is not None:
        add("CYP3A5 *3", f"{v} {_cp(v)} of *3", ["expresser (*1/*1)", "intermediate expresser (*1/*3)", "non-expresser (*3/*3)"][v],
            "tacrolimus", ["rs776746"], test_type="single variant (*3 only)", conseq=consequences_cyp3a5(v))
    v = _n(geno, "rs3745274")
    if v is not None:
        add("CYP2B6 *6", f"{v} {_cp(v)} of *6", ["no *6 detected", "intermediate metabolizer (*6 carrier)", "poor metabolizer (*6/*6)"][min(v, 2)],
            "efavirenz, methadone", ["rs3745274"], test_type="single variant (*6 only)", conseq=consequences_cyp2b6(v))
    v = _n(geno, "rs4148323")
    if v is not None:
        add("UGT1A1 *6", f"{v} {_cp(v)} of *6", ["no *6 detected", "reduced activity (*6 carrier)", "markedly reduced activity (*6/*6)"][min(v, 2)],
            "irinotecan, atazanavir", ["rs4148323"], test_type="single variant (*6 only; *28 repeat not tested)", conseq=consequences_ugt1a1(v))
    dp = [_n(geno, r) for r in ("rs3918290", "rs67376798", "rs55886062")]
    if None not in dp:
        add("DPYD", f"{sum(dp)} reduced-function alleles (*2A, D949V, *13 tested)",
            "no tested DPD-deficiency variant detected; residual risk from untested variants remains" if sum(dp) == 0
            else "reduced DPD activity: dose reduction is usually advised by guidelines",
            "fluorouracil, capecitabine", ["rs3918290", "rs67376798", "rs55886062"],
            test_type="limited screen (3 variants)", conseq=consequences_dpyd(sum(dp)))
    tp = [_n(geno, r) for r in ("rs1800462", "rs1800460", "rs1142345")]
    if None not in tp:
        add("TPMT", f"{sum(tp)} reduced-function alleles (*2, *3A/B/C tested)",
            (NEG if sum(tp) == 0 else "reduced TPMT activity"), "azathioprine, mercaptopurine, thioguanine",
            ["rs1800462", "rs1800460", "rs1142345"], test_type="limited screen (*2, *3B, *3C)", conseq=consequences_thiopurine(sum(tp), "TPMT"))
    v = _n(geno, "rs116855232")
    if v is not None:
        add("NUDT15 *3", f"{v} {_cp(v)}", (NEG if v == 0 else "reduced NUDT15 activity"),
            "thiopurines (azathioprine, mercaptopurine)", ["rs116855232"], test_type="limited screen (*3 only)", conseq=consequences_thiopurine(v, "NUDT15"))
    v = _n(geno, "rs1050828")
    if v is not None:
        add("G6PD (A- tag)", f"{v} {_cp(v)}", "A- variant (V68M) not detected; other G6PD deficiency variants are not assessed" if v == 0
            else "G6PD deficiency allele present",
            "rasburicase, primaquine, dapsone, nitrofurantoin (hemolysis risk if deficient)", ["rs1050828"],
            test_type="limited screen (one variant; low read depth is possible on X)", conseq=consequences_g6pd(v))
    v = _n(geno, "rs2395029")
    if v is not None:
        add("HLA-B*57:01 (tag only)", f"{v} {_cp(v)} of tag allele",
            "tag not detected (HLA typing is definitive)" if v == 0
            else "tag present: HLA-B*57:01 typing needed before relying on this (the tag is imperfect)",
            "abacavir (hypersensitivity if B*57:01 positive); flucloxacillin liver injury",
            ["rs2395029"], test_type="proxy SNP (not HLA typing)", conseq=consequences_hlab5701(v))
    return rows
