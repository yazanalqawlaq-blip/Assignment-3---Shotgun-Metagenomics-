#!/bin/bash
# running this script on login node

# =============================================================================
# 03_download_kraken2_db.sh
# Description: Download pre-built Kraken2 Standard-8 database
# =============================================================================

module load StdEnv/2023 gcc/12.3 kraken2/2.1.6

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
DBDIR=${PROJECT_DIR}/input_data/kraken2_db
mkdir -p ${DBDIR}
cd ${DBDIR}

# using the Standard-8 pre-built db (~8 GB)
# includes archaea, bacteria, viral, plasmid, human
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz

tar -xzf k2_standard_08_GB_20251015.tar.gz
rm k2_standard_08_GB_20251015.tar.gz

echo "done: $(date)"
ls -lh ${DBDIR}
