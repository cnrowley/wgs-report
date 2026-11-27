"""GeT-RM PGx benchmark: PharmCAT (via wgsreport.pgx_pharmcat) on 1000G 30x genomes vs CDC consensus diplotypes.

Run from a scratch directory (never inside the repo) that holds: truth/selected_samples.tsv (id + one column per gene, from the CDC
consolidated GeT-RM table), all60.vcf.gz (those samples, PharmCAT regions only, extracted from the NYGC 1000G 30x phased VCFs) and the
PharmCAT pipeline directory. See docs/VALIDATION.md for the commands. Public samples only."""
import json, subprocess, sys
from pathlib import Path
import pandas as pd
from wgsreport import pgx_pharmcat as ph

W = Path("run"); W.mkdir(exist_ok=True)
truth = pd.read_csv("truth/selected_samples.tsv", sep="\t").set_index("id")
rows = []
for sid in truth.index:
    out = W / sid
    summ = out / "summary.json"
    if not summ.exists():
        out.mkdir(exist_ok=True)
        vcf = out / f"{sid}.vcf.gz"
        subprocess.run(f"bcftools view -s {sid} all60.vcf.gz -Oz -o {vcf} && bcftools index -t {vcf}", shell=True, check=True)
        try:
            rep = ph.run(vcf, out, Path(sys.argv[1] if len(sys.argv) > 1 else "pharmcat"), None, sample=sid, absent_to_ref=True)
            summ.write_text(json.dumps(ph.load_summary(rep)))
        except Exception as e:
            print(sid, "FAILED", e, flush=True); continue
    s = json.loads(summ.read_text())
    calls = {g["gene"]: g for g in s["genes"]}
    for gene in ["CYP2C9", "CYP2C19", "TPMT", "CYP3A5", "SLCO1B1", "UGT1A1", "CYP2B6", "CYP4F2", "NUDT15", "CYP2D6"]:
        t = truth.at[sid, gene]
        if not isinstance(t, str):
            continue
        c = calls.get(gene)
        rows.append({"sample": sid, "gene": gene, "truth": t, "pharmcat": c["diplotype"] if c else "none", "ambiguous": bool(c and c["ambiguous"])})
    print(sid, "done", flush=True)
pd.DataFrame(rows).to_csv("bench_results.tsv", sep="\t", index=False)
