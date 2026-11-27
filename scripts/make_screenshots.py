#!/usr/bin/env python3
"""Make the documentation screenshots from a report built on a PUBLIC genome (never your own).

usage: make_screenshots.py RESULTS_DIR OUT_DIR
Cuts report.html into sections, renders each with headless Chrome/Chromium, trims the blank bottom. Needs google-chrome or chromium
and Pillow."""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image, ImageChops

res, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
chrome = next(c for c in ("google-chrome", "chromium", "chromium-browser") if shutil.which(c))
SHOTS = {  # name: (file, start marker, end marker)
    "01_ancestry": ("report.html", "<h2>1. Genetic ancestry", "<h2>2. Trait predictions"),
    "02_trait_scores": ("report.html", "<h2>2. Trait predictions", "<h3>Tier B: moderate"),
    "03_pharmacogenomics": ("report.html", "<h3>Pharmacogenomic summary", "<h3>Blood group"),
    "04_clinvar": ("report.html", "<h2>4. Rare diseases", "<h2>Inherited cancer"),
    "05_cancer": ("report.html", "<h2>Inherited cancer susceptibility", "<h2>4b. Carrier status"),
    "09_carrier": ("report.html", "<h2>4b. Carrier status", "<h2>Common clinical genetic tests"),
    "06_medication_handout": ("medication_report.html", "<h1>", "<h2>3. Drug lookup"),
    "08_clinical_tests": ("report.html", "<h2>Common clinical genetic tests", "<h2>Appendix"),
    "07_drug_lookup": ("medication_report.html", "<h2>3. Drug lookup", "<h2>4. Gene-by-gene"),
}


def shoot(html_path: Path, png: Path, width=1100):
    subprocess.run([chrome, "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars", f"--window-size={width},4200",
                    f"--screenshot={png}", f"file://{html_path}"], check=True, capture_output=True)
    im = Image.open(png).convert("RGB")
    bg = Image.new("RGB", im.size, im.getpixel((2, im.height - 2)))
    box = ImageChops.difference(im, bg).getbbox()
    im.crop((0, 0, im.width, min(im.height, (box[3] + 24) if box else im.height))).save(png)


for name, (fn, a, b) in SHOTS.items():
    src = (res / fn).read_text()
    i = src.find(a)
    j = src.find(b, i + 1) if i >= 0 else -1
    if i < 0 or j < 0:
        print("skip", name, "(marker not found)")
        continue
    head = src[:src.index("<body")] + "<body>"
    with tempfile.TemporaryDirectory() as t:
        p = Path(t) / "x.html"
        p.write_text(head + src[i:j] + "</body></html>")
        shoot(p, out / f"{name}.png")
    print("wrote", out / f"{name}.png")
