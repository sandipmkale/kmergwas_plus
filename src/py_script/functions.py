from glob import glob
import os
import sys
import subprocess
import shlex
import time

# Execution control (supports dry-run)
EXECUTE_COMMANDS = True

def set_dry_run(dry: bool):
    global EXECUTE_COMMANDS
    EXECUTE_COMMANDS = not dry

def get_file_type_in_dir(dir_name, type_suffix):
    fns = glob(f"{dir_name}/*.{type_suffix}")
    if len(fns) != 1:
        raise Exception(f"File count not = 1 {type_suffix} {len(fns)}\n{dir_name}")
    return fns[0]

def get_column(fn, index, sep="\t"):
    return [line.split(sep)[index] for line in open(fn, "r").read().split("\n")[:-1]]

def dir_exist(d):
    d = d[:-1] if len(d) > 0 and d[-1] == "/" else d
    return len(glob(d)) > 0

def create_dir(dir_path):
    if dir_exist(dir_path):
        raise Exception(f"directory already exist {dir_path}")
    os.makedirs(dir_path, exist_ok=True)
    return dir_path

def multiple_phenotypes_per_accession(fn):
    cmd = f"cat {shlex.quote(fn)} | tail -n +2 | cut -f1 | sort | uniq -c | awk '$1>1' | wc -l"
    try:
        out = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE).decode().strip()
        return int(out) > 0
    except Exception:
        return False

def run_and_log_command(cmd, logger):
    logger.writelines("RUN: " + cmd + "\n")
    if EXECUTE_COMMANDS:
        os.system(cmd)

def count_running_gemma(g_handles):
    return len([x for x in g_handles if x.poll() is None])

def run_gemma_cmd(cmd, maximal_to_run, g_handles, logger):
    while count_running_gemma(g_handles) >= maximal_to_run:
        time.sleep(30)
    if EXECUTE_COMMANDS:
        g_handles.append(subprocess.Popen(shlex.split(cmd)))
    logger.write(f"RUNG ({count_running_gemma(g_handles)}/{len(g_handles)}): {cmd}\n")

def copy_phenotypes(original_file, dest_file, average_script, logger):
    if multiple_phenotypes_per_accession(original_file):
        logger.writelines("Multiple phenotypes found per accessions -> Averaging\n")
        os.system(f"cat {original_file} | tail -n +2 | awk -f {average_script} > {dest_file}")
    else:
        logger.writelines("Unique phenotype per accession, copying phenotype data to directory\n")
        os.system(f"cp {original_file} {dest_file}")

def create_full_fam_file(new_fam_file, base_fam, fam_to_expand, n_phenotypes):
    orig_fam_acc = get_column(base_fam, 0, " ")
    cur_fam_acc = get_column(fam_to_expand, 0, "\t")[1:]
    cur_fam_pheno = [get_column(fam_to_expand, i+1, "\t")[1:] for i in range(n_phenotypes)]
    with open(new_fam_file, "w") as fout:
        for (i, acc_n) in enumerate(orig_fam_acc):
            fout.write(f"{acc_n} {acc_n} 0 0 0")
            if acc_n in cur_fam_acc:
                for p_ind in range(n_phenotypes):
                    fout.write(f" {cur_fam_pheno[p_ind][cur_fam_acc.index(acc_n)]}")
            else:
                for p_ind in range(n_phenotypes):
                    fout.write(" -9")
            fout.write("\n")

def get_file(fn):
    if len(glob(fn)) == 0:
        raise Exception(f"Couldn't find {fn}")
    return fn

def build_gemma_cmd(gemma_path, bed_prefix, kinship_fn, outdir, outname, maf, miss=0.5, covariates=None):
    bed = shlex.quote(bed_prefix + ".bed")
    fam = shlex.quote(bed_prefix + ".fam")
    bim = shlex.quote(bed_prefix + ".bim")
    cmd = (f"{shlex.quote(gemma_path)} -g {bed} -p {fam} -a {bim} "
           f"-lmm 2 -k {shlex.quote(kinship_fn)} "
           f"-outdir {shlex.quote(outdir)} -o {shlex.quote(outname)} -maf {maf} -miss {miss}")
    if covariates:
        cmd += f" -c {shlex.quote(covariates)}"
    return cmd

def create_new_relatedness_matrix(fn_in, fn_out, indices):
    """Subset and reorder a symmetric relatedness/kinship matrix (TSV) by 0-based indices."""
    # Read matrix
    with open(fn_in, "r") as f:
        rows = [line.rstrip().split() for line in f if line.strip()]
    # Convert to float
    try:
        mat = [[float(x) for x in row] for row in rows]
    except Exception as e:
        raise SystemExit(f"ERROR: Failed reading kinship matrix '{fn_in}': {e}")
    # Basic shape checks
    n = len(mat)
    for i, row in enumerate(mat):
        if len(row) != n:
            raise SystemExit(f"ERROR: Kinship matrix must be square; row {i} has length {len(row)} != {n}")
    # Subset rows and columns
    try:
        sub_rows = [mat[i] for i in indices]
        sub_mat = [[row[j] for j in indices] for row in sub_rows]
    except Exception as e:
        raise SystemExit(f"ERROR: Index error while subsetting kinship matrix: {e}")
    # Write out
    with open(fn_out, "w") as f:
        for row in sub_mat:
            f.write("\t".join(str(x) for x in row) + "\n")
