#!/usr/bin/env python3
"""
cleanup_intermediates_from_kmers_gwas.py

Safely remove intermediate files produced by kmers_gwas.py / kmers_gwas_mod.py,
mirroring the built-in cleanup logic from the original pipeline.

By default this script:
- Operates on subdirectories "kmers/" and/or "snps/" under --base-dir (if they exist)
- Removes temporary PLINK files created for permutations:
    kmers: pheno.*.P*.{bed,bim,fam,fam.orig}
    snps:  pheno.P*.{bed,bim,fam}
- Removes permutation outputs in each "output/" folder:
    rm output/P*
- Optionally gzips "output/phenotype_value.assoc.txt" to save space.

Use --dry-run to see what would be deleted or gzipped.
"""

import argparse
from pathlib import Path
import gzip
import shutil
import sys

def find_targets(base_dir: Path, targets_opt):
    present = []
    kmers_dir = base_dir / "kmers"
    snps_dir = base_dir / "snps"
    if (not targets_opt) or ("kmers" in targets_opt):
        if kmers_dir.exists():
            present.append(("kmers", kmers_dir))
    if (not targets_opt) or ("snps" in targets_opt):
        if snps_dir.exists():
            present.append(("snps", snps_dir))
    return present

def iter_globs(dirpath: Path, patterns):
    for pat in patterns:
        for p in dirpath.glob(pat):
            yield p

def gzip_if_exists(txt_path: Path, dry_run: bool, log):
    if not txt_path.exists():
        return
    gz_path = txt_path.with_suffix(txt_path.suffix + ".gz")
    if gz_path.exists():
        log.write(f"[skip] already gzipped: {gz_path}\n")
        return
    if dry_run:
        log.write(f"[dry-run] gzip {txt_path} -> {gz_path}\n")
        return
    # Stream gzip to avoid loading whole file
    with open(txt_path, "rb") as fin, gzip.open(gz_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    txt_path.unlink()
    log.write(f"[done] gzipped: {gz_path}\n")

def main():
    ap = argparse.ArgumentParser(description="Remove intermediate files from kmers_gwas outputs.")
    ap.add_argument("--base-dir", required=True, help="Directory containing 'kmers/' and/or 'snps/' subdirectories.")
    ap.add_argument("--targets", nargs="+", choices=["kmers", "snps"], help="Which subdirectories to clean. Default = all present.")
    ap.add_argument("--base-prefix", default="pheno", help="Base prefix used for PLINK perm files (default: pheno)")
    ap.add_argument("--gzip-main", action="store_true", help="Gzip output/phenotype_value.assoc.txt if present")
    ap.add_argument("--dry-run", action="store_true", help="List actions without making changes")
    args = ap.parse_args()

    base_dir = Path(args.base_dir).resolve()
    if not base_dir.exists():
        sys.exit(f"Base directory not found: {base_dir}")

    tasks = find_targets(base_dir, args.targets)
    if not tasks:
        sys.exit("No 'kmers' or 'snps' directories found or selected. Nothing to do.")

    log_path = base_dir / "cleanup_intermediates.log"
    with log_path.open("w") as log:
        log.write(f"Cleanup started in: {base_dir}\n")
        log.write(f"Targets: {', '.join([name for name,_ in tasks])}\n")
        log.write(f"Base prefix: {args.base_prefix}\n")
        log.write(f"Dry run: {args.dry_run}\n")

        for name, tdir in tasks:
            log.write(f"\n[{name}] in {tdir}\n")
            # Patterns mirroring kmers_gwas.py
            if name == "kmers":
                # pheno.*.P*.{bed,bim,fam,fam.orig}
                patterns = [
                    f"{args.base_prefix}.*.P*.bed",
                    f"{args.base_prefix}.*.P*.bim",
                    f"{args.base_prefix}.*.P*.fam",
                    f"{args.base_prefix}.*.P*.fam.orig",
                ]
                for p in iter_globs(tdir, patterns):
                    if args.dry_run:
                        log.write(f"[dry-run] rm {p}\n")
                    else:
                        try:
                            p.unlink()
                            log.write(f"[done] rm {p}\n")
                        except Exception as e:
                            log.write(f"[fail] rm {p}: {e}\n")
            elif name == "snps":
                # pheno.P*.{bed,bim,fam}
                patterns = [
                    f"{args.base_prefix}.P*.bed",
                    f"{args.base_prefix}.P*.bim",
                    f"{args.base_prefix}.P*.fam",
                ]
                for p in iter_globs(tdir, patterns):
                    if args.dry_run:
                        log.write(f"[dry-run] rm {p}\n")
                    else:
                        try:
                            p.unlink()
                            log.write(f"[done] rm {p}\n")
                        except Exception as e:
                            log.write(f"[fail] rm {p}: {e}\n")

            # Remove output/P*
            outdir = tdir / "output"
            if outdir.exists():
                for p in outdir.glob("P*"):
                    if args.dry_run:
                        log.write(f"[dry-run] rm {p}\n")
                    else:
                        try:
                            if p.is_dir():
                                # safety: do not recursively remove directories named P*
                                log.write(f"[skip] directory matched P* (not removing): {p}\n")
                            else:
                                p.unlink()
                                log.write(f"[done] rm {p}\n")
                        except Exception as e:
                            log.write(f"[fail] rm {p}: {e}\n")
                if args.gzip_main:
                    gzip_if_exists(outdir / "phenotype_value.assoc.txt", args.dry_run, log)

        log.write("\nCleanup completed.\n")

    print(f"Cleanup complete. Log saved to {log_path}")

if __name__ == "__main__":
    main()
