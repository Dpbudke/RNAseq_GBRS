# 3_CC_Implementation

Novel implementation of GBRS for Collaborative Cross (CC) and CCRIX mice. The standard GBRS tool was originally developed for Diversity Outbred (DO) mice; this directory contains the work required to extend it to CC and reciprocal intercross (CCRIX) populations, including generation of CC-specific transition probability matrices.

## Contents

### Scripts

| File | Description |
|------|-------------|
| `CC.trans.prob.Rmd` | Generates CC founder transition probability matrices used by GBRS for genotype reconstruction. Adapting this to CCRIX requires generating CCRIX-specific transprob files (in progress). |
| `cc_parser.py` | Helper script for parsing CC-format haplotype and concordance files. |

### Transition Probability Matrices

Required input for `1_GBRS_Pipeline/8_reconstruct.slurm`.

| File | Description |
|------|-------------|
| `tranprob.CC.G20.F.h5` | CC transition probability matrix — females, generation 20 (HDF5) |
| `tranprob.CC.G20.M.h5` | CC transition probability matrix — males, generation 20 (HDF5) |
| `tranprob.CC.G20.F.npz` | Same as above (NumPy compressed format) |
| `tranprob.CC.G20.M.npz` | Same as above (NumPy compressed format) |
| `RIX_CC001xCC019.F.h5` | CCRIX (CC001×CC019) transition probability — females |
| `RIX_CC019xCC001.F.h5` | CCRIX (CC019×CC001) transition probability — females |

### Reference Data (`Data/CC_indexes/`)

Haplotype concordance and strain index files used during CC implementation and validation.

| Directory | Contents |
|-----------|----------|
| `CC001/` | `concordance_results.csv`, `strain_index.csv` |
| `CC019/` | `concordance_results.csv`, `strain_index.csv` |
| `CC040/` | `concordance_results.csv`, `strain_index.csv` |
| (root) | `filtered_haplotypes_CC001/019/040.csv` |

## Usage

Open `CC.trans.prob.Rmd` in RStudio and knit to regenerate transition probability matrices:

```r
rmarkdown::render("CC.trans.prob.Rmd")
```

The output `.h5` files are then passed to `1_GBRS_Pipeline/8_reconstruct.slurm` via the `--tranprob-file` argument.

## Notes

- CCRIX transition probability generation is ongoing — the CC approach in `CC.trans.prob.Rmd` serves as the template to adapt for CCRIX crosses
