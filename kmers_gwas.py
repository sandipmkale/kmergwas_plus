#!/usr/bin/env python3
from __future__ import annotations
import os, sys, time, shlex
from pathlib import Path
from src.py_script.pipeline_parser import add_permutations_arg
paths: dict[str, str] = {}
ROOT = Path(__file__).resolve().parent
paths["associate_kmers"] = str(ROOT / "bin" / "associate_kmers")
paths["associate_snps"] = str(ROOT / "bin" / "associate_snps")
paths["gen_script"] = str(ROOT / "src" / "py_script" / "functions.py")
paths["parser_script"] = str(ROOT / "src" / "py_script" / "pipeline_parser.py")
paths["kinship_intersect_script"] = str(ROOT / "src" / "py_script" / "align_kinship_phenotype.py")
paths["Permute_phenotype"] = str(ROOT / "src" / "r_script" / "transform_and_permute_phenotypes.R")
paths["average_pheno_script"] = str(ROOT / "src" / "awk_scripts" / "average_phenotypes.awk")
paths["rscript_exec"] = "Rscript"
def _safe_exec_file(fp): 
    if fp and os.path.exists(fp):
        with open(fp,"rb") as fh: code = fh.read()
        exec(compile(code, fp, "exec"), globals(), globals())
_safe_exec_file(paths["gen_script"])
_safe_exec_file(paths["parser_script"])
def main(argv=None):
    import argparse
    local_parser = globals().get("parser")
    if local_parser is None:
        for fn_name in ("get_parser","build_parser","make_parser"):
            fn = globals().get(fn_name)
            if callable(fn):
                try:
                    local_parser = fn()
                    break
                except Exception:
                    pass
    if local_parser is None:
        p = argparse.ArgumentParser()
        p.add_argument("--outdir", default="out")
        p.add_argument("--gwas_software", choices=["gemma"], default="gemma")
        p.add_argument("--gemma_path"); p.add_argument("--blink_exec"); p.add_argument("--rscript_exec", default="Rscript")
        p.add_argument("--covariates")
        p.add_argument("--run_kmers", action="store_true"); p.add_argument("--run_one_step_snps", action="store_true"); p.add_argument("--run_two_steps_snps", action="store_true")
        p.add_argument("--fn_phenotype"); p.add_argument("--kmers_table"); p.add_argument("--snps_matrix")
        p.add_argument("--parallel", type=int, default=1); p.add_argument("--n_permutations", type=int, default=100)
        p.add_argument("--n_kmers", type=int, default=100000); p.add_argument("--n_snps", type=int, default=100000)
        p.add_argument("--kmers_len", type=int, default=31); p.add_argument("--maf", type=float, default=0.01); p.add_argument("--mac", type=int, default=5)
        p.add_argument("--use_kinship_from_kmers", action="store_true"); p.add_argument("--remove_intermediate", action="store_true"); p.add_argument("--min_data_points", type=int, default=10)
        local_parser = p
    args = local_parser.parse_args(argv)
    if getattr(args,"rscript_exec",None): paths["rscript_exec"]=args.rscript_exec
    globals()["args"] = args
    create_dir(args.outdir)  # noqa: F821
    paths["log_file"]=f"{args.outdir}/log_file"; f_log=open(paths["log_file"],"w",buffering=1); f_log.writelines(str(args)+"\n")
    base_created_files="pheno"
    paths["pheno_orig_fn"]=f"{args.outdir}/{base_created_files}.original_phenotypes"
    copy_phenotypes(get_file(args.fn_phenotype), paths["pheno_orig_fn"], paths["average_pheno_script"], f_log)  # noqa: F821
    if getattr(args,"snps_matrix",None):
        paths["snps_fam"]=args.snps_matrix+".fam"
    else:
        paths["snps_fam"]=f"{args.outdir}/{base_created_files}_for_accessions_order.fam"
        with open(paths["snps_fam"],"w") as fout:
            acc_order=[x for x in open(args.kmers_table+".names","r").read().split("\n") if x]
            for a in acc_order: fout.write(f"{a} {a} 0 0 0 -9\n")
    if args.use_kinship_from_kmers or (args.snps_matrix is None):
        paths["original_kinship_fn"]=args.kmers_table+".kinship"; f_log.write("Using kinship calculated on k-mers\n")
    else:
        paths["original_kinship_fn"]=args.snps_matrix+".kinship"; f_log.write("Using kinship calculated on SNPs\n")
    paths["pheno_intersected_fn"]=f"{args.outdir}/{base_created_files}.phenotypes"
    paths["kinship_fn"]=f"{args.outdir}/{base_created_files}.kinship"
    cur_cmd=(f"{shlex.quote(sys.executable)} {shlex.quote(paths['kinship_intersect_script'])} "
             f"--pheno {shlex.quote(paths['pheno_orig_fn'])} --fam_file {shlex.quote(paths['snps_fam'])} "
             f"--kinship_file {shlex.quote(paths['original_kinship_fn'])} --output_pheno {shlex.quote(paths['pheno_intersected_fn'])} "
             f"--output_kinship {shlex.quote(paths['kinship_fn'])} --DBs_list {shlex.quote(args.kmers_table+'.names')}")
    run_and_log_command(cur_cmd, f_log)  # noqa: F821
    paths["pheno_permuted_fn"]=f"{args.outdir}/{base_created_files}.phenotypes_and_permutations"
    paths["pheno_permuted_transformed_fn"]=f"{args.outdir}/{base_created_files}.phenotypes_permuted_transformed"
    paths["EMMA_perm_log_fn"]=f"{args.outdir}/EMMA_perm.log"
    paths["log_R_permute"]=f"{args.outdir}/phenotypes_transformation_permutation.log"
    cur_cmd=(f"{shlex.quote(paths['rscript_exec'])} {shlex.quote(paths['Permute_phenotype'])} "
             f"{shlex.quote(paths['pheno_intersected_fn'])} {shlex.quote(paths['kinship_fn'])} {int(args.n_permutations)} "
             f"{shlex.quote(paths['pheno_permuted_fn'])} {shlex.quote(paths['pheno_permuted_transformed_fn'])} "
             f"{shlex.quote(paths['EMMA_perm_log_fn'])} > {shlex.quote(paths['log_R_permute'])}")
    run_and_log_command(cur_cmd, f_log)  # noqa: F821
    pheno_file = paths.get("pheno_permuted_fn") if os.path.exists(paths.get("pheno_permuted_fn", "")) else paths["pheno_intersected_fn"]
    phenotypes_names = open(pheno_file, "r").read().split("\n")[0].split("\t")[1:]
    n_accession=len(open(paths["pheno_intersected_fn"]).read().split("\n")[1:-1])
    if n_accession < getattr(args,"min_data_points",10):
        f_log.write(f"Can't run with less than {args.min_data_points} data points; we have {n_accession}\n"); f_log.close(); Path(f"{args.outdir}/NOT_ENOUGH_DATA").touch(); return 0
    affective_maf=args.maf; maf_from_mac=float(args.mac)/float(n_accession); 
    if maf_from_mac>affective_maf: affective_maf=maf_from_mac
    print("affective MAF =", affective_maf)
    gemma_handles=[]
    if args.run_kmers:
        paths["kmers_associations_dir"]=f"{args.outdir}/kmers"; run_and_log_command(f"mkdir -p {paths['kmers_associations_dir']}", f_log)  # noqa: F821
        cur_cmd=(f"{shlex.quote(paths['associate_kmers'])} -p {shlex.quote(paths['pheno_permuted_transformed_fn'])} -b pheno -o {shlex.quote(paths['kmers_associations_dir'])} "
                 f"-n {int(args.n_kmers)} --parallel {int(args.parallel)} --kmers_table {shlex.quote(args.kmers_table)} --kmer_len {int(args.kmers_len)} --maf {affective_maf} --mac {int(args.mac)}")
        paths["log_kmers_associations"]=f"{args.outdir}/associate_kmers.log"; cur_cmd += f" 2> {shlex.quote(paths['log_kmers_associations'])}"
        run_and_log_command(cur_cmd, f_log)  # noqa: F821
        for (p_ind, p_name) in enumerate(phenotypes_names):
            fam_fn=get_file_type_in_dir(paths["kmers_associations_dir"], f"{p_name}.fam")  # noqa: F821
            run_and_log_command(f"mv {fam_fn} {fam_fn}.orig", f_log)  # noqa: F821
            base_name=fam_fn[:-4]
            cur_cmd = r"""cat %s | tail -n +2 | awk '{print $1 " " $1 " 0 0 0 " $%d}' > %s""" % (paths["pheno_permuted_fn"], p_ind+2, fam_fn)
            run_and_log_command(cur_cmd, f_log)  # noqa: F821
            covariates=getattr(args,"covariates",None);
            if not getattr(args,"gemma_path",None): raise SystemExit("ERROR: Provide --gemma_path /path/to/gemma")
            gemma_cmd=(f"{shlex.quote(args.gemma_path)} -bfile {shlex.quote(base_name)} -lmm 2 -k {shlex.quote(paths['kinship_fn'])} "
                           f"-outdir {shlex.quote(paths['kmers_associations_dir'] + '/output')} -o {shlex.quote(p_name)} -maf {affective_maf} -miss 0.5")
            if covariates: gemma_cmd += f" -c {shlex.quote(covariates)}"
            run_gemma_cmd(gemma_cmd, args.parallel, gemma_handles, f_log)  # noqa: F821
    while count_running_gemma(gemma_handles) != 0: time.sleep(1)  # noqa: F821
    f_log.close(); return 0
if __name__ == "__main__":
    sys.exit(main())