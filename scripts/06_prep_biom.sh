#!/bin/bash
#SBATCH --job-name=prep_biom
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=slurm/prep_biom_%j.out
#SBATCH --error=slurm/prep_biom_%j.err

# =============================================================================
# 06_prep_biom.sh
# Description: Convert Bracken reports to BIOM format for import into R
# Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
# Source: Human gut metagenome 
# =============================================================================

module load StdEnv/2023 python/3.11

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
BRDIR=${PROJECT_DIR}/output_files/bracken
OUTDIR=${PROJECT_DIR}/output_files/R_analysis
mkdir -p ${OUTDIR}

pip install kraken-biom --break-system-packages --quiet

# combined all bracken species reports into one BIOM file
# this is what i to import into phyloseq in R
kraken-biom \
    ${BRDIR}/SRR8146974_bracken_species.report \
    ${BRDIR}/SRR8146973_bracken_species.report \
    ${BRDIR}/SRR8146965_bracken_species.report \
    ${BRDIR}/SRR8146961_bracken_species.report \
    ${BRDIR}/SRR8146969_bracken_species.report \
    ${BRDIR}/SRR8146975_bracken_species.report \
    ${BRDIR}/SRR8146970_bracken_species.report \
    ${BRDIR}/SRR8146976_bracken_species.report \
    -o ${OUTDIR}/combined_bracken.biom

echo "=== done === $(date) ==="
ls -lh ${OUTDIR}/combined_bracken.biom
