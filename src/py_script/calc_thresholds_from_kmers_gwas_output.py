#!/usr/bin/env python3
"""
Calculate -log10(p) thresholds and filter significant associations from GWAS output.

Designed to post-process the output directory produced by kmers_gwas_mod.py (or compatible GWAS runs).
It reads permutation association result files, computes the empirical threshold(s) from the maxima
of -log10(p) across permutations, and filters the main association file by these thresholds.

Usage example:
  python3 calc_thresholds_from_kmers_gwas_output.py \
      --output-dir path/to/kmers_associations/output \
      --main-fn phenotype_value.assoc.txt \
      --alphas 0.05 0.10 \
      --pcol 10

Outputs:
  - threshold_5per, threshold_10per (one file per alpha) in the parent of --output-dir
  - pass_threshold_5per, pass_threshold_10per files with rows from main file exceeding threshold

Notes:
  - P-value column default is 10 (1-based) to match the original pipeline.
  - Supports both .assoc.txt and .assoc.txt.gz files.
  - Attempts to auto-detect a p-value header if present; falls back to --pcol if needed.
"""

from __future__ import annotations
import argparse, gzip, math, sys
from pathlib import Path
from typing import Iterable, List, Tuple, Optional

def open_maybe_gzip(fn: Path):
    if str(fn).endswith(".gz"):
        return gzip.open(fn, "rt")
    return open(fn, "rt")

def find_pval_index(header: str, default_idx_1based: int) -> int:
    """
    Given a header line, try to pick a p-value column 1-based index.
    Recognizes common header names: p, p_lrt, p_wald, p_score, pvalue, p.value, P, P_LRT, etc.
    If not found, returns default_idx_1based.
    """
    cols = header.strip().split()
    candidates = [c.lower() for c in cols]
    for i, c in enumerate(candidates, start=1):
        if c in {"p", "p_lrt", "pwald", "p_wald", "pvalue", "p.value", "pval", "p_score", "p_bonf", "p_fdr"}:
            return i
    # Sometimes GEMMA uses "p_lrt" exactly with tab delim; case-sensitive fallback
    for i, c in enumerate(cols, start=1):
        if c in {"p", "p_lrt", "p_wald", "pwald", "pvalue", "p.value", "pval", "p_score", "p_bonf", "p_fdr"}:
            return i
    return default_idx_1based

def parse_assoc_for_best_neglog10(fn: Path, default_pcol_1based: int, debug: bool=False) -> Optional[float]:
    """
    Return the maximum -log10(p) found in this association file.
    If file has a header, it will be used to detect p-value column.
    Returns None if no valid p-values are found.
    """
    best = None
    with open_maybe_gzip(fn) as f:
        first = f.readline()
        if not first:
            return None
        fields = first.strip().split()
        # if first line has non-numeric in the default p column, treat as header
        def try_float(x: str) -> bool:
            try:
                float(x)
                return True
            except Exception:
                return False

        pcol = default_pcol_1based
        if len(fields) >= default_pcol_1based and not try_float(fields[default_pcol_1based-1]):
            # header present; try to find p-value column by header name
            pcol = find_pval_index(first, default_pcol_1based)
            if debug:
                sys.stderr.write(f"[DEBUG] {fn.name}: header detected; using pcol={pcol}\n")
        else:
            # no header, process first line below
            val = float(fields[pcol-1])
            if val > 0:
                nl10 = -math.log10(val)
                best = nl10 if best is None else max(best, nl10)
            else:
                if debug:
                    sys.stderr.write(f"[DEBUG] {fn.name}: first line p<=0 ignored\n")

        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < pcol:
                if debug:
                    sys.stderr.write(f"[DEBUG] {fn.name}: line with insufficient columns skipped\n")
                continue
            try:
                p = float(parts[pcol-1])
                if p <= 0.0:
                    if debug:
                        sys.stderr.write(f"[DEBUG] {fn.name}: p<=0 encountered, skipping\n")
                    continue
                nl10 = -math.log10(p)
                best = nl10 if best is None else max(best, nl10)
            except Exception as e:
                if debug:
                    sys.stderr.write(f"[DEBUG] {fn.name}: parse error -> {e}\n")
                continue
    return best

def quantile_empirical(data: List[float], q: float) -> float:
    """
    Simple empirical quantile with linear interpolation between closest ranks.
    q in [0,1]. For threshold at alpha, we need the (1 - alpha) quantile of permutation maxima.
    """
    if not data:
        raise ValueError("No data provided for quantile computation.")
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

def write_text(fn: Path, text: str):
    fn.write_text(f"{text}\n")

def filter_main_by_threshold(main_fn: Path, out_fn: Path, threshold: float, default_pcol_1based: int):
    with open_maybe_gzip(main_fn) as fin, out_fn.open("w") as fout:
        header_checked = False
        pcol = default_pcol_1based
        for i, line in enumerate(fin):
            if not header_checked:
                fields = line.strip().split()
                # detect header
                try:
                    float(fields[default_pcol_1based-1])
                    has_header = False
                except Exception:
                    has_header = True
                if has_header:
                    pcol = find_pval_index(line, default_pcol_1based)
                    fout.write(line)  # keep header
                    header_checked = True
                    continue
                else:
                    header_checked = True
            parts = line.split()
            if len(parts) < pcol:
                if debug:
                    sys.stderr.write(f"[DEBUG] {fn.name}: line with insufficient columns skipped\n")
                continue
            try:
                p = float(parts[pcol-1])
            except Exception:
                continue
            if p > 0.0 and (-math.log10(p)) > threshold:
                fout.write(line)

def gather_assoc_files(output_dir: Path, main_basename: str) -> tuple[Path, list[Path]]:
    """
    Return main association file and list of permutation files from output_dir.
    Permutation files are all *.assoc.txt or *.assoc.txt.gz except the main file.
    """
    candidates = list(output_dir.glob("*.assoc.txt")) + list(output_dir.glob("*.assoc.txt.gz"))
    if not candidates:
        raise SystemExit(f"No association files (*.assoc.txt[.gz]) found in {output_dir}")
    # Helper to detect permutation files by name
    def is_perm(name: str) -> bool:
        n = name.lower()
        if "perm" in n or "permutation" in n or "shuffle" in n:
            return True
        import re
        return re.search(r"(?:^|[_\-\.])perm\d+(?:$|[_\-\.])", n) is not None
    # If a main file is explicitly given and not 'auto', honor it
    if main_basename and main_basename != "auto":
        main = None
        perms: list[Path] = []
        for fn in candidates:
            if fn.name == main_basename or fn.name == main_basename + ".gz":
                main = fn
            else:
                perms.append(fn)
        if main is None:
            raise SystemExit(f"Could not find main association file '{main_basename}' in {output_dir}")
        return main, perms
    # Auto-detect: choose non-permutation filename if unique; else fallback to latest modified
    non_perm = [c for c in candidates if not is_perm(c.name)]
    if len(non_perm) == 1:
        main = non_perm[0]
        perms = [c for c in candidates if c != main]
        return main, perms
    if len(non_perm) > 1:
        # Prefer a common default name if present
        preferred = ["phenotype_value.assoc.txt", "phenotype.assoc.txt", "trait.assoc.txt"]
        for p in preferred:
            for c in non_perm:
                if c.name == p or c.name == p + ".gz":
                    main = c
                    perms = [x for x in candidates if x != main]
                    return main, perms
        # Otherwise pick the most recently modified among non-perm files
        main = max(non_perm, key=lambda p: p.stat().st_mtime)
        perms = [x for x in candidates if x != main]
        return main, perms
    # If everything looks like a permutation by name, pick the most recent file as main
    main = max(candidates, key=lambda p: p.stat().st_mtime)
    perms = [x for x in candidates if x != main]
    return main, perms

def main():
    ap = argparse.ArgumentParser(description="Compute empirical -log10(p) thresholds from permutation outputs and filter main association file.")
    ap.add_argument("--output-dir", required=True, help="Path to the 'output' directory produced by kmers_gwas_mod.py")
    ap.add_argument("--main-fn", default="auto", help="Main association file name inside output-dir, or 'auto' to detect automatically (default: auto)")
    ap.add_argument("--alphas", type=float, nargs="+", default=[0.05, 0.10], help="Alpha levels for thresholds (default: 0.05 0.10)")
    ap.add_argument("--pcol", type=int, default=10, help="1-based index of p-value column if no header or auto-detect fails (default: 10)")
    ap.add_argument("--debug", action="store_true", help="Print diagnostics about file detection and p-value parsing")

    ap.add_argument("--out-parent", default=None, help="Directory to write threshold_* and pass_threshold_* files (default: parent of output-dir)")
    args = ap.parse_args()

    output_dir = Path(args.output_dir).resolve()
    if not output_dir.exists():
        raise SystemExit(f"Output directory not found: {output_dir}")
    out_parent = Path(args.out_parent).resolve() if args.out_parent else output_dir.parent

    log_fn = out_parent / "gwas_threshold_processing.log"
    with log_fn.open("w") as log:
        log.write("Starting threshold calculation process\n")
        log.write(f"Output directory: {output_dir}\n")
    main_fn, perm_files = gather_assoc_files(output_dir, args.main_fn)
    with log_fn.open("a") as log:
        log.write(f"Main association file detected: {main_fn.name}\n")
        log.write(f"Number of permutation files: {len(perm_files)}\n")
    if args.debug:
        sys.stderr.write(f"[DEBUG] Detected main: {main_fn}\n")
        sys.stderr.write(f"[DEBUG] Permutation files ({len(perm_files)}):\n" + "\n".join(["  - "+p.name for p in perm_files]) + "\n")

    # Collect permutation maxima of -log10(p)
    perm_maxima = []
    for fn in perm_files:
        best = parse_assoc_for_best_neglog10(fn, args.pcol, debug=args.debug)
        if best is not None:
            # Log best value
            with log_fn.open("a") as log:
                log.write(f"Permutation file {fn.name}: best -log10(p) = {best:.6f}\n")
            perm_maxima.append(best)
            if args.debug:
                sys.stderr.write(f"[DEBUG] {fn.name}: best -log10(p) = {best:.6f}\n")
        else:
            if args.debug:
                sys.stderr.write(f"[DEBUG] {fn.name}: no valid p-values found\n")

    if not perm_maxima:
        print("WARNING: No permutation files found or no valid p-values parsed. Thresholds cannot be computed.", file=sys.stderr)
        sys.exit(2)

    # Write best_pvals as -log10(p) values
    best_pvals_fn = out_parent / "best_pvals.txt"
    with best_pvals_fn.open("w") as bpf:
        for v in perm_maxima:
            bpf.write(f"{v:.6f}\n")
    with (out_parent / "gwas_threshold_processing.log").open("a") as log:
        log.write(f"Wrote best_pvals to {best_pvals_fn.name} (N={len(perm_maxima)})\n")

    # Compute thresholds at (1 - alpha) quantiles
    for alpha in args.alphas:
        q = 1 - alpha
        thr = quantile_empirical(perm_maxima, q)
        thr_fn = out_parent / f"threshold_{int(alpha*100)}per"
        write_text(thr_fn, f"{thr:.6f}")
        with log_fn.open("a") as log:
            log.write(f"Threshold for alpha={alpha}: {thr:.6f} written to {thr_fn.name}\n")
        # Filter main
        out_fn = out_parent / f"pass_threshold_{int(alpha*100)}per"
        filter_main_by_threshold(main_fn, out_fn, thr, args.pcol)
        with log_fn.open("a") as log:
            log.write(f"Filtered associations saved to {out_fn.name}\n")
        print(f"Wrote {thr_fn} and {out_fn}")
        if args.debug:
            sys.stderr.write(f"[DEBUG] Threshold for alpha={alpha}: {thr:.6f}\n")

if __name__ == "__main__":
    main()
