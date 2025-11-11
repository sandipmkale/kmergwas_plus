# Changelog

## [0.1.0] – 2025-11-03
### Added
- Initial refactored version of **kmersgwas-plus**.
- Modular directory structure (`src/py_script`, `r_script`, `awk_scripts`, `bin`).
- Unified runner `kmers_gwas.py` with direct GEMMA / BLINK / FarmCPUpp support.
- Covariate handling for GEMMA and FarmCPUpp.
- Complete README with usage examples and backend documentation.

### Fixed
- Removed obsolete `src/py/` path references.
- Corrected calls to AWK/R/Python helper scripts.
