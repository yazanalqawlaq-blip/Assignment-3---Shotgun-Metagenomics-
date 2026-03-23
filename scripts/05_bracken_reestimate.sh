#!/bin/bash
#SBATCH --job-name=bracken
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/bracken_%j.out
#SBATCH --error=slurm/bracken_%j.err

# =============================================================================
# 05_bracken_reestimate.sh
# Description: Bracken abundance re-estimation from Kraken2 reports
# Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
# Organism: Human gut metagenome 
# =============================================================================

module load StdEnv/2023 bracken/3.0

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
DBDIR=${PROJECT_DIR}/input_data/kraken2_db
K2DIR=${PROJECT_DIR}/output_files/kraken2
BRDIR=${PROJECT_DIR}/output_files/bracken
mkdir -p ${BRDIR}

SAMPLES="SRR8146974 SRR8146973 SRR8146965 SRR8146961 SRR8146969 SRR8146975 SRR8146970 SRR8146976"

# run bracken at species level for each sample
# this redistributes reads that kraken2 assigned at genus or higher
# back down to species, adjusting for genome size 
for SRR in ${SAMPLES}; do
    echo "--- ${SRR} --- $(date) ---"

    bracken \
        -d ${DBDIR} \
        -i ${K2DIR}/${SRR}.k2report \
        -o ${BRDIR}/${SRR}.bracken \
        -w ${BRDIR}/${SRR}_bracken_species.report \
        -r 150 \
        -l S \
        -t 10

    # -r 150 because our reads are ~150 bp 
    # -l S for species level
    # -t 10 filters out anything with fewer than 10 reads

    echo "${SRR} done"
done

echo "=== bracken complete === $(date) ==="
ls -lh ${BRDIR}/*.bracken
