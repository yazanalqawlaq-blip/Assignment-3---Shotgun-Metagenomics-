#!/bin/bash
# no SLURM headers here since compute nodes don't have internet
# have to run this directly on the login node

# =============================================================================
# 01_download_data.sh
# Description: Download raw shotgun metagenomics reads from NCBI SRA
# Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
# Organism: Human gut metagenome (Italian cohort)
# =============================================================================

module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
RAWDIR=${PROJECT_DIR}/input_data/raw_reads
mkdir -p ${RAWDIR}
cd ${RAWDIR}

# sra-toolkit fills up home with cache if you don't redirect it
mkdir -p /scratch/yazanalq/sra_cache
echo '/repository/user/main/public/root = "/scratch/yazanalq/sra_cache"' > ~/.ncbi/user-settings.mkfg

# 4 vegan + 4 omnivore, all paired-end NextSeq 500
# chose these for similar depth (~9-11 GB) so comparisons aren't biased
# vegans are from a mix of Turin/Bari/Parma, omnivores mostly Turin
VEGAN="SRR8146974 SRR8146973 SRR8146965 SRR8146961"
OMNIVORE="SRR8146969 SRR8146975 SRR8146970 SRR8146976"
ALL="${VEGAN} ${OMNIVORE}"

for SRR in ${ALL}; do
    echo "--- ${SRR} --- $(date) ---"

    # skip if already downloaded
    if [ -f "${SRR}_1.fastq.gz" ] && [ -f "${SRR}_2.fastq.gz" ]; then
        echo "already exists, skipping"
        continue
    fi

    # prefetch first then convert, way more reliable for big files
    prefetch ${SRR} --max-size 20G -O .
    fasterq-dump ${SRR} --split-files --threads 4 --progress

    # compress immediately, uncompressed fastqs are massive
    gzip ${SRR}_1.fastq
    gzip ${SRR}_2.fastq

    # clean up the .sra file so we don't waste space
    rm -rf ${SRR}

    echo "${SRR} done"
done

echo "=== all downloads complete === $(date) ==="
ls -lh ${RAWDIR}/*.fastq.gz
