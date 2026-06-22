#!/bin/bash
#SBATCH --job-name=kcs_sweep
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=1-12%10
#SBATCH --output=logs_kraken/kcs_%A_%a.out
#SBATCH --error=logs_kraken/kcs_%A_%a.err

# Kraken2 Confidence-Score sweep -> karioCaS-ready *.mpa files.
#
# Key change vs the original wrapper: output files are named SAMPLE_CS<value>.mpa
# (was *_conf_*), which is what karioCaS's importer parses. The decimal form is
# unambiguous: _CS0.0 -> 0%, _CS0.1 -> 10%, _CS0.5 -> 50%, _CS1.0 -> 100%.
# Files are written straight into <PROJECT>/000_mpa_original/.

module load kraken2/2.17.1

DB="/opt/share/databases/Kraken2_Bracken_RefSeq_Genomes_Complete"

# Paired reads to classify (one *_R1*.fastq.gz per sample).
READS_DIR="$HOME/data.thiago.parente/karioCaS_validation/reads"
# karioCaS project root; the .mpa files land in $PROJECT/000_mpa_original/
PROJECT="$HOME/data.thiago.parente/karioCaS_validation/kariocas_project"
OUT_DIR="${PROJECT}/000_mpa_original"

mkdir -p "${OUT_DIR}" logs_kraken

R1_FILES=("${READS_DIR}"/*_R1.fastq.gz)
INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))
FASTQ_R1="${R1_FILES[$INDEX]}"
FASTQ_R2="${FASTQ_R1/_R1.fastq.gz/_R2.fastq.gz}"

# Sample name: basename minus _R1.fastq.gz. MUST be free of underscores if you
# want group overlays (group = name with trailing digits stripped). For the mock
# use names like MOCKeven01, MOCKstag01, ...
SAMPLE_NAME=$(basename "${FASTQ_R1}" _R1.fastq.gz)

echo "Sample: ${SAMPLE_NAME}"
echo "  R1: ${FASTQ_R1}"
echo "  R2: ${FASTQ_R2}"

for CONFIDENCE in $(seq 0.0 0.1 1.0); do
    CONF=$(printf "%.1f" "${CONFIDENCE}")
    REPORT_FILE="${OUT_DIR}/${SAMPLE_NAME}_CS${CONF}.mpa"   # e.g. MOCKeven01_CS0.5.mpa
    echo "  -> Kraken2 paired, confidence ${CONF} -> $(basename "${REPORT_FILE}")"
    kraken2 \
        --db "${DB}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory-mapping \
        --paired \
        --confidence "${CONF}" \
        --use-mpa-style \
        --report "${REPORT_FILE}" \
        --output /dev/null \
        "${FASTQ_R1}" "${FASTQ_R2}"
done

echo "Done: ${SAMPLE_NAME}"
