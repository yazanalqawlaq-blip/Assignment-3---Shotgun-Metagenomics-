#!/bin/bash
#SBATCH --job-name=qc_fastp
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm/qc_fastp_%j.out
#SBATCH --error=slurm/qc_fastp_%j.err

# =============================================================================
# 02_qc_fastp.sh
# Description: QC and adapter trimming on raw paired-end shotgun reads
# Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
# Organism: Human gut metagenome (Italian cohort)
# =============================================================================

module load StdEnv/2023 fastp/1.0.1

PROJECT_DIR=/scratch/yazanalq/Assignment-3---Shotgun-Metagenomics-
RAWDIR=${PROJECT_DIR}/input_data/raw_reads
QCDIR=${PROJECT_DIR}/output_files/qc
mkdir -p ${QCDIR}

SAMPLES="SRR8146974 SRR8146973 SRR8146965 SRR8146961 SRR8146969 SRR8146975 SRR8146970 SRR8146976"

for SRR in ${SAMPLES}; do
    echo "--- running fastp on ${SRR} --- $(date) ---"

    fastp \
        -i ${RAWDIR}/${SRR}_1.fastq.gz \
        -I ${RAWDIR}/${SRR}_2.fastq.gz \
        -o ${QCDIR}/${SRR}_trimmed_1.fastq.gz \
        -O ${QCDIR}/${SRR}_trimmed_2.fastq.gz \
        --html ${QCDIR}/${SRR}_fastp_report.html \
        --json ${QCDIR}/${SRR}_fastp_report.json \
        --thread 8 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --detect_adapter_for_pe \
        --correction \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20

    # phred 20 = 99% base call accuracy, pretty standard cutoff
    # length_required 50 keeps reads long enough to be useful for kmer classification
    # cut_front/cut_tail does sliding window trimming from both ends
    # --correction fixes mismatched bases in overlapping paired-end regions
    # --detect_adapter_for_pe auto detects and removes illumina adapters

    echo "${SRR} done"
done

echo "=== fastp complete === $(date) ==="

# quick summary of what we kept vs filtered
for SRR in ${SAMPLES}; do
    echo "--- ${SRR} ---"
    grep -A 2 '"filtering_result"' ${QCDIR}/${SRR}_fastp_report.json | head -5
done
