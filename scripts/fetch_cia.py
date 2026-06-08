#!/usr/bin/env python3
"""Download the HITRAN collision-induced-absorption (CIA) files used for Titan's
thermal-IR opacity, then build the compact band-averaged table.

Sources (HITRAN CIA, Karman et al. 2019 and updates; hitran.org/cia):
  N2-N2  : main/N2-N2_2021.cia    (0-450 cm^-1, 70-400 K)
  N2-CH4 : main/N2-CH4_2024.cia   (0-800 cm^-1, 70-400 K)
  N2-H2  : main/N2-H2_2024.cia    (0-1886 cm^-1)
  CH4-CH4: alternate/CH4-CH4_2011.cia

The raw .cia files (~13 MB) are written to data/cia/raw/ (gitignored); the
small band-averaged table data/cia/titan_cia_bands.npz is committed.

Usage:  python3 scripts/fetch_cia.py
"""

import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RAW = ROOT / "data" / "cia" / "raw"
BASE = "https://hitran.org/data/CIA/"
FILES = {
    "N2-N2": "main/N2-N2_2021.cia",
    "N2-CH4": "main/N2-CH4_2024.cia",
    "N2-H2": "main/N2-H2_2024.cia",
    "CH4-CH4": "alternate/CH4-CH4_2011.cia",
}
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"


def fetch():
    RAW.mkdir(parents=True, exist_ok=True)
    for pair, rel in FILES.items():
        out = RAW / Path(rel).name
        if out.exists():
            print(f"{pair}: have {out.name}")
            continue
        print(f"{pair}: downloading {rel} ...")
        req = urllib.request.Request(BASE + rel, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            out.write_bytes(r.read())
        print(f"  -> {out} ({out.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    fetch()
    print("\nBuilding band-averaged table ...")
    from rt.cia import build_band_table
    build_band_table()
    print("done.")
