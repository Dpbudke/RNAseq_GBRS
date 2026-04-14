# RNAseq_GBRS

Genome Reconstruction by RNA-Seq (GBRS) pipeline for allele-specific expression (ASE) quantification in Collaborative Cross (CC), CC-RIX, and Diversity Outbred (DO) mice. GBRS leverages the known founder haplotype structure of CC mice to reconstruct local genotypes from RNA-seq reads and quantify founder strain-specific expression.

## Overview

GBRS uses RNA-seq reads to simultaneously perform:
1. **Genotype reconstruction** — infer founder haplotype probabilities at each genomic locus from read alignment patterns
2. **Allele-specific expression** — quantify expression separately for each of the 8 CC founder alleles

This approach is based on [Choi et al. 2021 (GBRS)](https://github.com/churchill-lab/gbrs) and the [EMASE](https://github.com/churchill-lab/emase) framework.

## Repository Structure

```
RNAseq_GBRS/
├── 1_GBRS_Pipeline/      # HPC SLURM scripts: alignment → quantification → genotype reconstruction
└── 2_GBRS_Analysis/      # R Markdown: transition probability modeling, ASE analysis
```

## Workflow

```
Preprocessed FASTQs (from RNAseq_STAR)
    │
    ▼
1_GBRS_Pipeline
    ├── Bowtie2 alignment to CC diploid reference
    ├── BAM sort/index
    ├── BAM → EMASE format
    ├── Multi-way alignment merge
    ├── EMASE quantification (allele-specific counts)
    ├── Genotype reconstruction (GBRS reconstruct)
    ├── Allele-specific expression quantification
    ├── Count combination and confirmation
    ├── Genotype interpolation
    └── Genotype export
    │
    ▼
2_GBRS_Analysis
    ├── Transition probability matrix generation (CC founders)
    └── ASE visualization and interpretation
```

## Dependencies

| Stage | Tools |
|-------|-------|
| HPC pipeline | GBRS, Bowtie2, SAMtools, EMASE, SLURM |
| R analysis | R, GBRS R package, ape, ggplot2 |

## Data Files (in `2_GBRS_Analysis/`)

| File | Description |
|------|-------------|
| `tranprob.CC.G20.F.h5` / `.M.h5` | CC founder transition probability matrices (female/male), generation 20 |
| `tranprob.CC.G20.F.npz` / `.M.npz` | Same as above in NumPy compressed format |
| `RIX_CC001xCC019.F.h5` | CCRIX transition probability file |
| `Data/CC_indexes/` | Haplotype concordance and strain index files for CC001, CC019, CC040 |

## Related Repository

Standard STAR alignment and DE analysis pipeline: [RNAseq_STAR](https://github.com/Dpbudke/RNAseq_STAR)

## Author

Dawson Budke — dpbudke@ucdavis.edu
