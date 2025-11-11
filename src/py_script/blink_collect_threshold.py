#!/usr/bin/env python3
import os
import sys
import math
import gzip
import argparse
from pathlib import Path
from typing import Optional, List, Tuple

def open_maybe_gzip(fn: Path):
    return gzip.open(fn, "rt") if str(fn).endswith(".gz") else open(fn, "rt")

def find_pcol_from_header(header: str, default_idx_1based: int = 1) -> int:
    cols = header.strip().split()
    low = [c.lower() for c in cols]
    for i, c in enumerate(low, start=1):
        if c in {"p_value", "p", "pvalue", "p.val", "p_lrt", "pwald", "p_wald"}:
            return i
    return default_idx_1based

def best_stats_from_file(fn: Path, default_pcol_1based: int = 1) -> Optional[Tuple[float, float]]:
    best_neg = None
    best_p = None
    with open_maybe_gzip(fn) as f:
        first = f.readline()
        if not first:
            return None
        cols = first.strip().split()
        try:
            float(cols[default_pcol_1based-1])
            header = False
            pcol = default_pcol_1based
        except Exception:
            header = True
            pcol = find_pcol_from_header(first, default_pcol_1based)
        if not header:
            try:
                p = float(cols[pcol-1])
                if p > 0:
                    neg = -math.log10(p)
                    best_neg = neg if best_neg is None else max(best_neg, neg)
                    best_p = p if best_p is None or p < best_p else best_p
            except Exception:
                pass
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < pcol:
                continue
            try:
                p = float(parts[pcol-1])
                if p <= 0:
                    continue
                neg = -math.log10(p)
                best_neg = neg if best_neg is None else max(best_neg, neg)
                best_p = p if best_p is None or p < best_p else best_p
            except Exception:
                continue

    if best_neg is None:
        return None
    return best_p, best_neg

def empirical_quantile(data: List[float], q: float) -> float:
    xs = sorted(data)
    if q <= 0: return xs[0]
    if q >= 1: return xs[-1]
    pos = (len(xs)-1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    frac = pos - lo
    return xs[lo] * (1-frac) + xs[hi] * frac

def collect_main_files(input_dir: Path) -> List[Path]:
    return sorted([p for p in input_dir.glob("*.txt") if p.name.lower().endswith("_pheno_gwas_result.txt")])

def main():
    ap = argparse.ArgumentParser(description="Compute thresholds from multiple BLINKC result files and record best p-values.")
    ap.add_argument("--input-dir", required=True, help="Directory containing BLINKC result files")
    ap.add_argument("--alphas", type=float, nargs="+", default=[0.05, 0.10], help="Alpha levels for thresholds")
    args = ap.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        print(f"ERROR: directory not found: {input_dir}", file=sys.stderr)
        sys.exit(2)

    main_files = collect_main_files(input_dir)
    if not main_files:
        print("No BLINKC outputs found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(main_files)} BLINKC result files.")

    best_rows = []  # (file, best_p, best_neglog10)

    for fn in main_files:
        stats = best_stats_from_file(fn, default_pcol_1based=1)
        if stats is None:
            print(f"[WARN] No valid p-values in {fn.name}")
            continue
        best_p, best_neg = stats
        best_rows.append((fn, best_p, best_neg))
        print(f"{fn.name}: best p={best_p:.6g}, best -log10(p)={best_neg:.6f}")

    # ✅ Write summary of best p-values
    summary_file = input_dir / "blinkc_best_pval"
    with summary_file.open("w") as out:
        out.write("run_name\tbest_pval\tbest_neglog10p\n")
        for fn, best_p, best_neg in best_rows:
            out.write(f"{fn.stem}\t{best_p:.6g}\t{best_neg:.6f}\n")

    print(f"\nWrote best p-value summary → {summary_file}")

    # ✅ Compute thresholds
    maxima = [neg for (_, _, neg) in best_rows]
    for alpha in args.alphas:
        thr = empirical_quantile(maxima, 1 - alpha)
        thr_file = input_dir / f"threshold_{int(alpha*100)}per"
        with open(thr_file, "w") as out:
            out.write(f"{thr:.6f}\n")
        print(f"Wrote {thr_file} = {thr:.6f}")

if __name__ == "__main__":
    main()

