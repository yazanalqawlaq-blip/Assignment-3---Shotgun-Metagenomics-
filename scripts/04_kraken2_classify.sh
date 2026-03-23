#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=slurm/kraken2_%j.out
#SBATCH --error=slurm/kraken2_%j.err

# =============================================================================
# 04_kraken2_classify.sh
# Description: Taxonomic classification of QC'd reads using Kraken2
# Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
# Organism: Human gut metagenome 
# =============================================================================

module load StdEnv/2023 gcc/12.3 kraken2/2.1.6

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
DBDIR=${PROJECT_DIR}/input_data/kraken2_db
QCDIR=${PROJECT_DIR}/output_files/qc
K2DIR=${PROJECT_DIR}/output_files/kraken2
mkdir -p ${K2DIR}

SAMPLES="SRR8146974 SRR8146973 SRR8146965 SRR8146961 SRR8146969 SRR8146975 SRR8146970 SRR8146976"

for SRR in ${SAMPLES}; do
    echo "--- ${SRR} --- $(date) ---"

    kraken2 \
        --db ${DBDIR} \
        --paired \
        --gzip-compressed \
        --confidence 0.15 \
        --threads 8 \
        --output ${K2DIR}/${SRR}.kraken2 \
        --report ${K2DIR}/${SRR}.k2report \
        ${QCDIR}/${SRR}_trimmed_1.fastq.gz \
        ${QCDIR}/${SRR}_trimmed_2.fastq.gz

    echo "${SRR} done"
done

# confidence 0.15 as recommended in lecture to reduce false positives
# my data is paired-end and gzipped so i needed to add flags to get the script to run

echo "=== kraken2 complete === $(date) ==="
ls -lh ${K2DIR}/*.k2report
