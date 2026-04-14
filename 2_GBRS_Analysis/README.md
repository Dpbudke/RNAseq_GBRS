# 2_GBRS_Analysis

R Markdown scripts for GBRS transition probability modeling and allele-specific expression (ASE) visualization in CC and CCRIX mice.

## Scripts

### `CC.trans.prob.Rmd`
Generates and validates transition probability matrices for CC founder strains used by GBRS for genotype reconstruction.

- **Output:** `tranprob.CC.G20.F.h5`, `tranprob.CC.G20.M.h5`, `.npz` equivalents
- **Requires:** CC founder haplotype data, R GBRS package

### `GBRS_ASE.Rmd`
Visualizes and analyzes allele-specific expression results from the GBRS pipeline. Includes strain-level expression summaries and ASE plots.

### `cc_parser.py`
Helper Python script for parsing CC-format haplotype/concordance files.

## Data Files

| File | Description |
|------|-------------|
| `tranprob.CC.G20.F.h5` | CC transition probability matrix — females, generation 20 |
| `tranprob.CC.G20.M.h5` | CC transition probability matrix — males, generation 20 |
| `tranprob.CC.G20.F.npz` | Same as above (NumPy format) |
| `tranprob.CC.G20.M.npz` | Same as above (NumPy format) |
| `RIX_CC001xCC019.F.h5` | CCRIX (CC001×CC019) transition probability — females |
| `RIX_CC019xCC001.F.h5` | CCRIX (CC019×CC001) transition probability — females |
| `Data/CC_indexes/` | Per-strain haplotype concordance and strain index files (CC001, CC019, CC040) |

## Usage

Open any `.Rmd` in RStudio and knit, or:

```r
rmarkdown::render("CC.trans.prob.Rmd")
rmarkdown::render("GBRS_ASE.Rmd")
```

## Notes

- The `.h5` transition probability files generated here are required as input by `1_GBRS_Pipeline/8_reconstruct.slurm`
- CCRIX transition probability generation is ongoing — see `CC.trans.prob.Rmd` for the CC approach to adapt for CCRIX
