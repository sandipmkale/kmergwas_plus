#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pathlib import Path

def find_fams(root: Path):
    return sorted(root.rglob("*.fam"))

def make_pheno(fam_path: Path) -> Path:
    """
    Convert .fam → phenotype text file (same directory, same base name).
    Output format:
        taxa    pheno
        IID     phenotype_value
    """
    rows = []
    with fam_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            parts = ln.split()
            if len(parts) >= 6:
                iid = parts[1]
                pheno = parts[5]
                rows.append((iid, pheno))

    if not rows:
        raise RuntimeError(f"No valid phenotype entries found in: {fam_path}")

    # pheno.5.P5.fam → pheno.5.P5.txt
    pheno_txt = fam_path.with_suffix(".txt")

    with pheno_txt.open("w", encoding="utf-8") as out:
        out.write("taxa\tpheno\n")
        for iid, phval in rows:
            out.write(f"{iid}\t{phval}\n")

    return pheno_txt

def run_cmd(cmd):
    print("\n[RUN]", " ".join(map(str, cmd)))
    return subprocess.run(cmd).returncode

def main():
    ap = argparse.ArgumentParser(description="Generate phenotype files from .fam and run BLINKC.")
    ap.add_argument("--blinkc", required=True, help="Path to BLINKC executable")
    ap.add_argument("--input-dir", required=True, help="Directory containing .fam files")
    ap.add_argument("--output-dir", required=True, help="Directory for BLINKC output")
    args = ap.parse_args()

    blinkc = Path(args.blinkc).resolve()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    fams = find_fams(input_dir)
    if not fams:
        print(f"ERROR: No .fam files found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(fams)} .fam files in {input_dir}")

    for fam in fams:
        pheno_txt = make_pheno(fam)

        base = fam.stem  # e.g., pheno.5.P5
        pheno_base = str(pheno_txt.with_suffix(""))  # remove .txt extension

        out_prefix = output_dir / base

        cmd = [
            str(blinkc),
            "--gwas",
            "--file", pheno_base,
            "--plink",
            "--out", str(out_prefix)
        ]

        rc = run_cmd(cmd)
        if rc != 0:
            print(f"BLINKC ERROR on {fam}", file=sys.stderr)

    print("\nDone.\n")

if __name__ == "__main__":
    main()

