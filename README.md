# kmersGWAS+
A K-mer–based Genome-Wide Association Analysis Framework Integrating GEMMA and BLINK-C

## Origin and Attribution

This repository is based on the kmersGWAS pipeline originally developed by Yotam Voichek in the laboratory of Detlef Weigel (Max Planck Institute for Biology Tübingen). The original pipeline demonstrated the power of k-mer–based GWAS in genetically diverse genomes and provided reference-independent genotype representations.

Original repository:
https://github.com/voichek/kmersGWAS

Original publication:
Voichek, Y., & Weigel, D. (2020). Identifying genetic variants underlying phenotypic variation without complete genomes. Nature Biotechnology, 38(5), 560–564.
https://doi.org/10.1038/s41587-020-0441-4

Please cite the original authors when using or adapting this work.

## Overview

kmersGWAS+ extends the original workflow by adding support for BLINK-C association modeling, phenotype batching, optional covariates, empirical threshold estimation, and post-run storage optimization.

The pipeline performs genome-wide association studies directly on k-mer presence/absence matrices, providing reference-independent genotyping support suitable for organisms with high genomic structural variation.

## Key Features
- Support for GEMMA (LMM) and BLINK-C association testing
- Reference-free k-mer genotype matrices
- Optional kinship correction
- Empirical p-value threshold estimation from permutations
- Multi-trait batch processing
- Output cleanup utilities

## Dependencies

| Software | Purpose |
|---------|---------|
| Python ≥ 3.8 | Pipeline logic |
| GEMMA | Mixed-model association testing |
| BLINK-C | Bayesian model selection–based GWAS |
| PLINK ≥ 1.9 | Phenotype and sample structure utilities |
| NumPy, SciPy, Pandas | Core computation libraries |

## Inputs

1. k-mer matrix (binary or count-based)
2. Phenotypes (FAM or tabular format)
3. Optional covariates
4. Optional permutation output directories

## Script Descriptions

| Script | Location | Purpose |
|-------|----------|---------|
| kmers_gwas.py | project root | Main GWAS pipeline orchestrating k-mer matrix handling, kinship estimation (optional), and GEMMA/BLINK-C execution. |
| functions.py | src/py_script/ | Shared helper routines for matrix handling, phenotype alignment, and association computation. |
| pipeline_parser.py | src/py_script/ | Defines command-line interface settings and run parameters. |
| calc_thresholds_from_kmers_gwas_output.py | src/py_script/ | Computes empirical significance thresholds from permutation results. |
| cleanup_intermediates_from_kmers_gwas.py | src/py_script/ | Removes large temporary output files to manage storage use. |
| run_gwas_using_blinkc.py | src/py_script/ | Generates phenotype tables and executes BLINK-C across traits. |
| blink_collect_threshold.py | src/py_script/ | Aggregates minimal p-values from BLINK-C permutations to compute empirical thresholds. |

## Usage Examples

### BLINK-C GWAS
```bash
python src/py_script/run_gwas_using_blinkc.py     --blinkc /path/to/blinkc     --input-dir path/to/plink_fam_directory     --output-dir results/blinkc_output
```

### Threshold Collection
```bash
python src/py_script/blink_collect_threshold.py     --input-dir results/blinkc_output     --output results/blinkc_thresholds.txt
```

### Full Pipeline (GEMMA)
```bash
python kmers_gwas.py     --kmer-matrix data/kmers_matrix.tsv     --phenotype data/phenotypes.fam     --model gemma     --output results/gemma_output     --compute-kinship
```

## Output

- Per-trait association statistics
- Optional kinship matrix
- Empirical significance thresholds (if permutation workflow used)
- Optional reduced summary results

## Citation

If publishing results, please cite:

Voichek, Y., & Weigel, D. (2020). Identifying genetic variants underlying phenotypic variation without complete genomes. Nature Biotechnology, 38(5), 560–564.

## License

Specify license details here (e.g., MIT, GPL-3.0, CC-BY 4.0).
