"""`wgsreport analyze`: score one personal VCF and assemble the report inputs."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from . import disclaimer, provenance, build, ancestry, clinvar, measured as measured_mod, bloodgroup, pgs, pgx, rare, refpanel, score, snps, variant_scores, carrier
from .registry import score_ids, for_ancestry, load_registry
from .report import write_report
from . import medication_report
from .utils import log, require_tools, run

FILE_FLOOR_COVERAGE = 0.80   # warn if less of a score's variants could be matched
PC_ADJUST = 4                # number of PCs regressed out within the reference population


def _load_copynumber(path, psychiatric: bool):
    """Output of `copy-number`; `--exclude psychiatric` keeps the CNV screen but drops its neuropsychiatric estimates."""
    if not path:
        return None
    cn = json.loads(Path(path).read_text())
    if not psychiatric and cn.get("cnv_screen"):
        cn["cnv_screen"]["psychiatric_estimates"] = False
    return cn


def detect_sex(vcf: Path, requested: str) -> tuple[str, str]:
    if requested != "auto":
        return requested, "user-specified"
    hdr = run(["bcftools", "view", "-h", str(vcf)]).stdout
    for line in hdr.splitlines():
        if line.startswith("##referenceSexKaryotype="):
            k = line.split("=")[1].strip()
            return ("male" if k == "XY" else "female" if k == "XX" else "unknown"), "VCF header karyotype"
    # heterozygosity on non-PAR X (GRCh38: X:2.78-155.7 Mb)
    for chrom in ("chrX", "X"):
        try:
            out = run(["bcftools", "view", "-H", "-v", "snps", "-r", f"{chrom}:3000000-150000000", str(vcf)]).stdout
        except RuntimeError:
            continue
        het = tot = 0
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 10:
                continue
            alleles = f[9].split(":")[0].replace("|", "/").split("/")   # GT subfield only; haploid calls (males on non-PAR X) have one allele
            if "." in alleles:
                continue
            tot += 1
            het += len(alleles) == 2 and alleles[0] != alleles[1]
        if tot:
            return ("female" if het / tot > 0.25 else "male"), f"chrX heterozygosity {het / tot:.2f}"
    return "unknown", "undetermined"


def vcf_summary(vcf: Path) -> dict:
    out = run(["bcftools", "stats", str(vcf)]).stdout
    d = {}
    for l in out.splitlines():
        if l.startswith("SN\t"):
            _, _, k, v = l.split("\t")
            d[k.strip(": ")] = int(v)
        elif l.startswith("TSTV\t"):
            d["ts/tv"] = float(l.split("\t")[4])
    return d


def _z(raw: float, ref: pd.Series, pcs_ref: pd.DataFrame | None, user_pcs: np.ndarray | None):
    """Return (z_raw, z_pc_adjusted) of `raw` relative to reference scores `ref`."""
    z_raw = (raw - ref.mean()) / ref.std(ddof=1)
    if pcs_ref is None:
        return z_raw, z_raw
    X = np.column_stack([np.ones(len(ref)), pcs_ref.loc[ref.index].to_numpy()[:, :PC_ADJUST]])
    beta, *_ = np.linalg.lstsq(X, ref.to_numpy(), rcond=None)
    resid = ref.to_numpy() - X @ beta
    sd = resid.std(ddof=X.shape[1])
    pred = np.concatenate([[1.0], user_pcs[:PC_ADJUST]]) @ beta
    return z_raw, (raw - pred) / sd


def callability_check(dos: np.ndarray, assumed: np.ndarray, p: np.ndarray) -> dict:
    """Compare observed genotype-class fractions with Hardy-Weinberg expectations.

    A variant-only VCF cannot prove that an absent site is homozygous reference. If many absent sites were
    really non-callable, the observed hom-ref fraction would exceed the expectation from reference allele
    frequencies (and het/hom-alt fractions would fall short). Matching fractions are reassuring, not proof.
    """
    ok = ~np.isnan(p) & ~np.isnan(dos)
    p, d = p[ok], dos[ok]
    n = len(d)
    exp = {"homref": float(((1 - p) ** 2).mean()), "het": float((2 * p * (1 - p)).mean()), "homalt": float((p ** 2).mean())}
    obs = {"homref": float((d == 0).mean()), "het": float((d == 1).mean()), "homalt": float((d == 2).mean())}
    bins = {}
    for lo, hi in ((0, .05), (.05, .2), (.2, .5), (.5, .8), (.8, 1.0001)):
        m = (p >= lo) & (p < hi)
        if m.sum():
            bins[f"{lo:.2f}-{min(hi, 1):.2f}"] = {"n": int(m.sum()), "obs_homref": float((d[m] == 0).mean()),
                                                  "exp_homref": float(((1 - p[m]) ** 2).mean())}
    return {"n": n, **{f"obs_{k}": v for k, v in obs.items()}, **{f"exp_{k}": v for k, v in exp.items()}, "bins": bins}


MAX_CENTROID_SD = 3.0   # beyond this the genome is off the reference cloud: standardised scores may be biased (Ashkenazi genomes, ~3.4 SD, showed inflated height z)


def analyze(args) -> None:
    log(disclaimer.TEXT.replace("`", ""))
    require_tools("bcftools", "plink2")
    names = run(["bcftools", "query", "-l", str(args.vcf)]).stdout.split()
    if len(names) != 1:      # genotypes of several samples would be concatenated and silently misread
        raise SystemExit(f"The VCF has {len(names)} samples; wgs-report needs exactly one. Extract yours first: "
                         "bcftools view -s SAMPLE -o one.vcf.gz -Oz in.vcf.gz")
    if not getattr(args, "skip_build_check", False):
        err = build.check_build(run(["bcftools", "view", "-h", str(args.vcf)]).stdout)
        if err:
            raise SystemExit(err)
    work, out = Path(args.work_dir), Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    reg = load_registry()
    traits = [t for t in reg["traits"] if not t.get("sensitive") or t["sensitive"] in args.sensitive]

    log("loading reference data")
    weights = pd.read_pickle(work / "weights.pkl")
    pvar = pd.read_pickle(work / "pvar.pkl")
    af_all = pd.read_pickle(work / "af_all.pkl")
    panel_scores = pd.read_csv(work / "panel_scores.tsv", sep="\t", index_col=0)
    panel_pcs = pd.read_csv(work / "panel_pcs.tsv", sep="\t").set_index("IID")
    meta = pd.read_csv(work / "panel" / "samples.tsv", sep="\t")
    files = json.load(open(work / "pca_files.json"))
    scale = np.load(work / "pc_scale.npy")
    sstats = json.load(open(work / "score_stats.json"))["stats"]

    traits = [t for t in traits if t["pgs"] in panel_scores.columns]
    sex, sex_src = detect_sex(args.vcf, args.sex)
    log(f"sex: {sex} ({sex_src})")

    # ---- positions we need from the user's VCF --------------------------------------------
    snp_panel = snps.load_panel("apoe" in args.sensitive)
    plp = None
    if not args.no_clinvar:
        cv = args.clinvar or clinvar.fetch(work / "clinvar.vcf.gz")
        log("loading ClinVar pathogenic/likely-pathogenic records")
        min_stars = getattr(args, 'clinvar_min_stars', 2)
        plp_all = clinvar.load_plp(cv, min_stars=1)
        plp = plp_all[plp_all["stars"] >= min_stars]
        # single-star P/LP records in cancer-predisposition genes are shown separately, as lower review status, not dropped silently
        cg = clinvar.cancer_genes()
        plp_low = plp_all[(plp_all["stars"] < min_stars) & plp_all["gene"].isin(cg)]
        carrier_low_set = plp_all[(plp_all["stars"] < min_stars) & plp_all["gene"].isin(carrier.recessive_genes())]
        log(f"  {len(plp):,} P/LP records (>={min_stars} review stars); {len(plp_low):,} lower-review records in cancer genes")
    positions = set(zip(pvar["chrom"], pvar["pos"].astype(int))) | snps.panel_positions(snp_panel) | bloodgroup.POSITIONS
    variant_sets = {}
    for t in traits:
        if t.get("variant_set"):
            vdf, _ = pgs.read_score(pgs.fetch_score(t["variant_set"], work / "pgs_cache"))
            variant_sets[t["id"]] = vdf
            positions |= variant_scores.positions(vdf)
    if plp is not None:
        positions |= set(zip(plp["chrom"], plp["pos"]))
        positions |= set(zip(plp_low["chrom"], plp_low["pos"]))
        positions |= set(zip(carrier_low_set["chrom"], carrier_low_set["pos"]))
    log(f"reading user VCF at {len(positions):,} positions")
    sites = score.read_user_sites(args.vcf, positions)

    # ---- genotypes -> scores ----------------------------------------------------------------
    log("computing genotype dosages")
    dos, dstats, assumed = score.user_dosages(sites, pvar)
    if not getattr(args, "skip_build_check", False):
        err = build.check_build("", dstats)
        if err:
            raise SystemExit(err)
    w_aligned = weights.reindex(pvar["id"]).fillna(0.0)
    mean_alt = af_all.reindex(pvar["id"]).to_numpy()
    raw = score.user_scores(dos, w_aligned, mean_alt)
    callability = callability_check(dos, assumed, refpanel.read_afreq(work / "panel", "EUR").reindex(pvar["id"]).to_numpy())
    log(f"callability check: observed hom-ref {callability['obs_homref']:.3f} vs expected {callability['exp_homref']:.3f}")

    # ---- ancestry -----------------------------------------------------------------------------
    log("ancestry projection")
    user_pcs = ancestry.project_user(dos, pvar, files) * scale
    pcs_cols = [c for c in panel_pcs.columns if c.startswith("PC")]
    anc = ancestry.classify(user_pcs, panel_pcs.reset_index(), meta)
    group = anc["superpop"]
    ref_ids = meta.loc[meta.super_pop == group, "IID"]
    ref_ids = [i for i in ref_ids if i in panel_scores.index]
    ancestry_ok = anc["superpop_fractions"].get(group, 0) >= 0.9 and anc["distance_from_centroid_sd"] < MAX_CENTROID_SD
    log(f"ancestry: {group} ({anc['superpop_fractions']}) centroid distance {anc['distance_from_centroid_sd']:.1f} SD")

    # ---- per-trait results ------------------------------------------------------------------------
    results, dist_data = [], {}
    avail = set(panel_scores.columns)
    for t in traits:
        t = for_ancestry(t, group, ancestry_ok, avail)
        pid = t["pgs"]
        st = sstats[pid]
        if t.get("sex") and t["sex"] != sex:
            results.append({**t, "skipped": f"applies to {t['sex']}s only"})
            continue
        ref = panel_scores.loc[ref_ids, pid]
        z_raw, z = _z(raw[pid], ref, panel_pcs.loc[ref_ids, pcs_cols], user_pcs)
        res = {**t, "skipped": None, "z_raw": float(z_raw), "z": float(z),
               "percentile": float(100 * norm.cdf(z)), "coverage": st["coverage"],
               "n_scored": st["n_scored"], "n_variants": st["n_variants"],
               "low_coverage": st["coverage"] < FILE_FLOOR_COVERAGE, "n_ref": len(ref_ids)}
        if t.get("crosscheck"):
            cc = t["crosscheck"]
            _, zc = _z(raw[cc], panel_scores.loc[ref_ids, cc], panel_pcs.loc[ref_ids, pcs_cols], user_pcs)
            res["crosscheck_z"] = float(zc)
            res["crosscheck_pgs"] = cc
        if t["kind"] == "quantitative":
            r = np.sqrt(t["r2"])
            half = 1.96 * np.sqrt(1 - t["r2"])
            res.update(expected_sd=float(r * z), pi_sd=(float(r * z - half), float(r * z + half)),
                       expected_percentile=float(100 * norm.cdf(r * z)))
            if t["tier"] != "C" and t.get("unit") and sex in (t.get("mean") or {}):
                mu, sd = t["mean"][sex], t["sd"][sex]
                if t["id"] == "height" and args.height_mean:
                    mu = args.height_mean
                if t["id"] == "height" and args.height_sd:
                    sd = args.height_sd
                res.update(abs_unit=t["unit"], abs_mean=mu, abs_sd=sd,
                           abs_expected=float(mu + sd * r * z),
                           abs_interval=(float(mu + sd * (r * z - half)), float(mu + sd * (r * z + half))))
        else:
            res["relative_odds"] = float(t["or_per_sd"] ** z)
        if t.get("variant_set"):
            vs = variant_scores.score_variant_set(variant_sets[t["id"]], sites)
            res["variant_set"] = {"pgs": t["variant_set"], "unit": t.get("variant_set_unit", ""), **vs}
            vid = t["variant_set"]
            if vid in panel_scores.columns and vid in raw.index:       # where this total falls among reference individuals
                refv = panel_scores.loc[ref_ids, vid]
                res["variant_set"]["reference"] = {"group": group, "n": int(len(refv)),
                                                   "fraction_at_or_above": float((refv >= raw[vid] - 1e-6).mean()),
                                                   "fraction_any_carried": float((refv > 1e-6).mean())}
        tc = score.top_contributors(dos, w_aligned[pid], pvar, mean_alt, n=6)
        tc["contribution_sd"] = tc["contribution"] / ref.std(ddof=1)
        res["top_variants"] = tc.drop(columns=["rsid_hint"]).to_dict("records")
        results.append(res)
        dist_data[t["id"]] = (ref.to_numpy(), float(raw[pid]), float(ref.mean()), float(ref.std(ddof=1)), z)

    # ---- SNP panel, ClinVar, QC -------------------------------------------------------------------------
    cram_json = json.load(open(args.cram_json)) if args.cram_json else None
    geno = snps.genotype(snp_panel, sites, cram_json["panel"] if cram_json else None)
    apoe = snps.decode_apoe(geno) if "apoe" in args.sensitive else None
    cv_hits = clinvar.screen(plp, sites) if plp is not None else pd.DataFrame()
    cancer_low = clinvar.screen(plp_low, sites) if plp is not None else pd.DataFrame()
    carrier_low = clinvar.screen(carrier_low_set, sites) if plp is not None else pd.DataFrame()
    log(f"ClinVar P/LP hits: {len(cv_hits)}")
    summary = vcf_summary(args.vcf)
    pgx_rows = pgx.summarize(geno, args.cyp2d6)
    lof = pd.DataFrame()
    if args.lof_scan:
        if not args.fasta:
            raise SystemExit("--lof-scan needs --fasta")
        log("scanning ACMG, carrier-panel, cancer-predisposition and AF genes for rare predicted loss-of-function variants")
        (work / "rare").mkdir(exist_ok=True)
        af_genes = rare.load_af_genes()
        lof = rare.lof_scan(args.vcf, args.fasta, rare.prepare_gff(work / "rare"),
                            (clinvar.acmg_genes() | set(af_genes) | set(clinvar.carrier_genes()) | set(clinvar.cancer_genes())) - rare.NO_LOF_MECHANISM, args.kg_dir)
        lof = rare.drop_clinvar_matches(lof, cv_hits)
        if len(lof):
            lof["list"] = [", ".join(x for x in (("ACMG secondary findings" if g in clinvar.acmg_genes() else ""),
                                                 (f"AF gene set (OR {af_genes[g]})" if g in af_genes else ""),
                                                 (f"ACMG carrier panel ({clinvar.carrier_genes()[g].rsplit('(', 1)[-1].rstrip(')')})" if g in clinvar.carrier_genes() else ""),
                                                 (f"ClinGen cancer gene ({clinvar.cancer_genes()[g]})" if g in clinvar.cancer_genes() else "")) if x)
                           for g in lof["gene"]]
        log(f"  {len(lof)} candidate variants")

    meas_rows = measured_mod.compare(results, measured_mod.load_wide(args.measured), args.treated_from) if args.measured else []
    out_json = {
        "vcf": args.vcf.name, "sex": sex, "sex_source": sex_src, "ancestry": anc, "ancestry_ok": ancestry_ok,
        "reference_group": group, "n_reference": len(ref_ids), "user_pcs": user_pcs.tolist(),
        "dosage_stats": dstats, "callability": callability, "vcf_summary": summary, "apoe": apoe,
        "pgx": pgx_rows, "carrier_low_review": carrier_low.drop(columns=["acmg_sf_gene", "carrier_panel_gene", "cancer_gene"], errors="ignore").to_dict("records"),
        "cancer_low_review": cancer_low.drop(columns=["acmg_sf_gene", "carrier_panel_gene"], errors="ignore").to_dict("records"),
        "pharmcat": json.load(open(args.pharmcat_json)) if getattr(args, "pharmcat_json", None) else None,
        "repeats": json.load(open(args.repeats_json)) if getattr(args, "repeats_json", None) else None,
        "y_haplogroup": json.load(open(args.y_haplogroup_json)) if getattr(args, "y_haplogroup_json", None) else None,
        "qc": json.load(open(args.qc_json)) if getattr(args, "qc_json", None) else None,
        "hla": json.load(open(args.hla_json)) if getattr(args, "hla_json", None) else None, "abo": bloodgroup.infer_abo(sites), "mito": (cram_json or {}).get("mito"), "lof_scan_run": bool(args.lof_scan),
        "lof_candidates": lof.to_dict("records"),
        "admixture": json.load(open(args.admixture_json)) if args.admixture_json else None,
        "copynumber": _load_copynumber(args.copynumber_json, "psychiatric" in args.sensitive),
        "mt_haplogroup": args.mt_haplogroup, "measured": meas_rows,
        "traits": [{k: v for k, v in r.items() if k != "top_variants"} for r in results],
        "disclaimer": disclaimer.TEXT,
        "provenance": provenance.collect(sorted(score_ids(reg)), work / "pgs_cache", work,
                                         {"build_check_skipped": bool(getattr(args, "skip_build_check", False))}),
    }
    json.dump(out_json, open(out / "results.json", "w"), indent=1, default=float)
    geno.to_csv(out / "snp_panel_genotypes.tsv", sep="\t", index=False)
    cv_hits.to_csv(out / "clinvar_pathogenic_hits.tsv", sep="\t", index=False)
    lof.to_csv(out / "rare_lof_candidates.tsv", sep="\t", index=False)
    write_report(out, out_json, results, dist_data, geno, cv_hits, panel_pcs, meta, user_pcs, reg)
    log(f"report written to {out / 'report.html'}")
    log(f"printable medication report: {medication_report.write(out, out_json)}")
