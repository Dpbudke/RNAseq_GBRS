# 1_GBRS_Pipeline

SLURM-based HPC scripts for the full GBRS allele-specific expression pipeline in CC mice. Scripts are designed to run sequentially as SLURM jobs on an HPC cluster.

## Pipeline Steps

| Script | Step | Description |
|--------|------|-------------|
| `1_gbrs_alignment.slurm` | 1 | Bowtie2 alignment to CC diploid transcriptome |
| `2_bowtie_stats.slurm` | 2a | Summarize Bowtie2 alignment statistics |
| `2_combine_txt_files.sh` | 2b | Combine alignment stat text files |
| `3_bam_sort_index.slurm` | 3 | Sort and index BAM files (SAMtools) |
| `4_bam2emase.slurm` | 4 | Convert BAM to EMASE alignment format |
| `5_common_align_merge.slurm.py` | 5 | Merge multi-way EMASE alignments |
| `6_compress_emase.slurm` | 6 | Compress EMASE alignment files |
| `7_quantify.slurm` | 7 | EMASE quantification (allele-specific counts) |
| `8_reconstruct.slurm` | 8 | GBRS genotype reconstruction |
| `9_quantify_ASE.slurm` | 9 | Allele-specific expression quantification |
| `10_combine_counts.py` | 10 | Combine per-sample count files |
| `11_confirm_reads2counts.py` | 11 | Validate count matrix integrity |
| `12_interpolate.slurm` | 12 | Interpolate genotype probabilities |
| `13_plot_genome.slurm` | 13 | Plot genome-wide genotype reconstruction |
| `14_genotype_export.slurm` | 14 | Export reconstructed genotypes |

## Usage

Run scripts in numerical order. Each script requires `samples.txt` (one sample name per line) and appropriate reference files (CC diploid transcriptome, transition probability matrices).

```bash
mkdir -p slurmout
sbatch 1_gbrs_alignment.slurm
# Wait for completion, then:
sbatch 3_bam_sort_index.slurm
# etc.
```

## Notes

- Scripts were developed for CC pup ASE analysis; extending to CCRIX requires generation of CCRIX-specific transition probability files (see `2_GBRS_Analysis/`)
- Transition probability matrices for CC mice (generation 20) are provided in `2_GBRS_Analysis/`
