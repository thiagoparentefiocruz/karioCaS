#!/bin/bash
#SBATCH --job-name=camisim
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=logs_camisim/camisim_%j.out
#SBATCH --error=logs_camisim/camisim_%j.err
set -euo pipefail
mkdir -p logs_camisim

# Adjust to your CAMISIM checkout / conda env.
CAMISIM=/path/to/CAMISIM
# source activate camisim   # or: conda activate camisim

python "${CAMISIM}/metagenomesimulation.py" validation/camisim_config.ini

# CAMISIM writes paired FASTQs per sample under camisim_out/<sample>/reads/.
# Stage them as <SAMPLE>_R1.fastq.gz / <SAMPLE>_R2.fastq.gz for the Kraken sweep:
mkdir -p reads
while read -r S; do
  r1=$(find camisim_out -path "*${S}*" -name '*1.fq.gz' -o -path "*${S}*" -name '*1.fastq.gz' | head -1)
  r2=$(find camisim_out -path "*${S}*" -name '*2.fq.gz' -o -path "*${S}*" -name '*2.fastq.gz' | head -1)
  [ -n "${r1:-}" ] && ln -sf "$(readlink -f "$r1")" "reads/${S}_R1.fastq.gz"
  [ -n "${r2:-}" ] && ln -sf "$(readlink -f "$r2")" "reads/${S}_R2.fastq.gz"
done < camisim_in/samples.txt
echo "Reads staged in reads/  ->  next: kraken_cs_sweep.sbatch"

# ---------------------------------------------------------------------------
# FALLBACK if CAMISIM setup stalls: InSilicoSeq (pip install InSilicoSeq) gives
# the same paired Illumina reads from our genomes + abundances, per sample:
#   iss generate --genomes <concat all .fna> --abundance_file <per-sample> \
#       --model hiseq --n_reads 5M --cpus 24 --output reads/<SAMPLE>
# (InSilicoSeq abundance_file = genome_id <tab> fraction; reuse our distributions
#  after mapping G### ids to the concatenated FASTA headers.)
