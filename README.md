# RNAseq_GBRS

Genome Reconstruction by RNA-Seq (GBRS) pipeline for allele-specific expression (ASE) quantification in Collaborative Cross (CC), CC-RIX, and Diversity Outbred (DO) mice. GBRS leverages the known founder haplotype structure of CC mice to reconstruct local genotypes from RNA-seq reads and quantify founder strain-specific expression.

A key contribution of this repository is the novel implementation of GBRS for CC and CCRIX mice (`3_CC_Implementation/`), including generation of CC-specific transition probability matrices required for genotype reconstruction.

## Overview

GBRS uses RNA-seq reads to simultaneously perform:
1. **Genotype reconstruction** — infer founder haplotype probabilities at each genomic locus from read alignment patterns
2. **Allele-specific expression** — quantify expression separately for each of the 8 CC founder alleles

This approach is based on [Choi et al. 2021 (GBRS)](https://github.com/churchill-lab/gbrs) and the [EMASE](https://github.com/churchill-lab/emase) framework.

## Repository Structure

```
RNAseq_GBRS/
├── 1_GBRS_Pipeline/       # HPC SLURM scripts: alignment → quantification → genotype reconstruction
├── 2_GBRS_Analysis/       # R Markdown: ASE visualization and interpretation
└── 3_CC_Implementation/   # Novel CC/CCRIX transition probability generation and data files
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
2_GBRS_Analysis             ← uses transition probability files from 3_CC_Implementation/
    └── ASE visualization and interpretation
```

## Dependencies

| Stage | Tools |
|-------|-------|
| HPC pipeline | GBRS, Bowtie2, SAMtools, EMASE, SLURM |
| R analysis | R, GBRS R package, ape, ggplot2 |

## Reference Data

Reference files required for EMASE and GBRS (diploid transcriptome, founder genomes, etc.) are archived on Zenodo:

> Lloyd, Michael. (2025). *Reference Data for EMASE and GBRS*. Zenodo. https://doi.org/10.5281/zenodo.15091310

Dawson Budke is a contributor to this dataset.

## Related Repository

Standard STAR alignment and DE analysis pipeline: [RNAseq_STAR](https://github.com/Dpbudke/RNAseq_STAR)

## Author

Dawson Budke — dpbudke@ucdavis.edu
