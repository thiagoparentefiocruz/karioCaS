#!/bin/bash
# Download a multi-domain RefSeq genome set for the mock community via the NCBI
# `datasets` CLI (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
#   conda install -c conda-forge ncbi-datasets-cli   # or `module load` on HPC
#
# Output:
#   genomes/<safe_name>.fna        one FASTA per genome
#   genome_manifest.tsv            file  taxid  organism  domain   (used by build_community.py)
set -euo pipefail

OUT=genomes
mkdir -p "${OUT}" tmp_dl
: > genome_manifest.tsv
printf "file\ttaxid\torganism\tdomain\n" >> genome_manifest.tsv

# domain <TAB> species name  (all present in RefSeq with reference genomes)
SPECIES=$(cat <<'EOF'
Bacteria	Escherichia coli
Bacteria	Staphylococcus aureus
Bacteria	Pseudomonas aeruginosa
Bacteria	Bacillus subtilis
Bacteria	Klebsiella pneumoniae
Bacteria	Salmonella enterica
Bacteria	Listeria monocytogenes
Bacteria	Enterococcus faecalis
Bacteria	Streptococcus pneumoniae
Bacteria	Clostridioides difficile
Bacteria	Mycobacterium tuberculosis
Bacteria	Helicobacter pylori
Bacteria	Campylobacter jejuni
Bacteria	Vibrio cholerae
Bacteria	Neisseria meningitidis
Bacteria	Haemophilus influenzae
Bacteria	Acinetobacter baumannii
Bacteria	Bacteroides fragilis
Bacteria	Bifidobacterium longum
Bacteria	Corynebacterium glutamicum
Bacteria	Caulobacter vibrioides
Bacteria	Deinococcus radiodurans
Bacteria	Thermus thermophilus
Bacteria	Treponema pallidum
Bacteria	Chlamydia trachomatis
Bacteria	Borreliella burgdorferi
Bacteria	Lactobacillus delbrueckii
Bacteria	Streptomyces coelicolor
Bacteria	Fusobacterium nucleatum
Bacteria	Porphyromonas gingivalis
Archaea	Methanocaldococcus jannaschii
Archaea	Halobacterium salinarum
Archaea	Saccharolobus solfataricus
Archaea	Methanosarcina barkeri
Archaea	Pyrococcus furiosus
Eukaryota	Saccharomyces cerevisiae
Eukaryota	Candida albicans
Eukaryota	Cryptococcus neoformans
Eukaryota	Aspergillus fumigatus
Eukaryota	Plasmodium falciparum
Viruses	Escherichia phage Lambda
Viruses	Enterobacteria phage T4
Viruses	Human betaherpesvirus 5
Viruses	Influenza A virus
Viruses	Vaccinia virus
EOF
)

while IFS=$'\t' read -r DOMAIN NAME; do
    [ -z "${NAME:-}" ] && continue
    safe=$(echo "${NAME}" | tr ' ' '_' | tr -cd 'A-Za-z0-9_')
    echo ">>> ${DOMAIN}: ${NAME}"
    rm -rf tmp_dl/pkg && mkdir -p tmp_dl/pkg
    # Reference RefSeq genome + sequence report (for taxid).
    datasets download genome taxon "${NAME}" \
        --reference --assembly-source RefSeq \
        --include genome,seq-report \
        --filename "tmp_dl/${safe}.zip" || { echo "  [skip] download failed"; continue; }
    unzip -oq "tmp_dl/${safe}.zip" -d tmp_dl/pkg
    fna=$(find tmp_dl/pkg -name '*.fna' | head -1)
    [ -z "${fna}" ] && { echo "  [skip] no FASTA"; continue; }
    cp "${fna}" "${OUT}/${safe}.fna"
    # taxid from the data report
    taxid=$(dataformat tsv genome \
        --package "tmp_dl/${safe}.zip" \
        --fields organism-tax-id 2>/dev/null | tail -n +2 | head -1)
    taxid=${taxid:-NA}
    printf "%s\t%s\t%s\t%s\n" "${OUT}/${safe}.fna" "${taxid}" "${NAME}" "${DOMAIN}" >> genome_manifest.tsv
done <<< "${SPECIES}"

# Host contaminant (not part of truth). Use one chromosome to stay light, or
# swap for the full GRCh38 if you prefer a heavier host load.
echo ">>> Host: Homo sapiens (GRCh38, chr21 as contaminant)"
datasets download genome accession GCF_000001405.40 \
    --chromosomes 21 --include genome --filename tmp_dl/human.zip || true
if [ -f tmp_dl/human.zip ]; then
    unzip -oq tmp_dl/human.zip -d tmp_dl/human
    hfna=$(find tmp_dl/human -name '*.fna' | head -1)
    [ -n "${hfna}" ] && cp "${hfna}" "${OUT}/Homo_sapiens_chr21.fna" \
        && printf "%s\t9606\tHomo sapiens\tHost\n" "${OUT}/Homo_sapiens_chr21.fna" >> genome_manifest.tsv
fi

echo "Done. $(($(wc -l < genome_manifest.tsv) - 1)) genomes -> ${OUT}/ ; manifest: genome_manifest.tsv"
