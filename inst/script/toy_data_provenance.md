# Provenance of the example data in `inst/extdata`

`inst/extdata/your_project_name/000_mpa_original/` contains six Kraken2
reports in MPA format for a single sample, `SAMPLE01`, classified at six
Kraken2 confidence scores (CS):

| File                 | Kraken2 `--confidence` |
|----------------------|------------------------|
| `SAMPLE01_CS00.mpa`  | 0.0                    |
| `SAMPLE01_CS02.mpa`  | 0.2                    |
| `SAMPLE01_CS04.mpa`  | 0.4                    |
| `SAMPLE01_CS06.mpa`  | 0.6                    |
| `SAMPLE01_CS08.mpa`  | 0.8                    |
| `SAMPLE01_CS09.mpa`  | 0.9                    |

## Source

The reads come from one of the author's own shotgun metagenomic samples
(unpublished; not deposited in a public archive at the time of submission).

## How the files were produced

1. The sequencing reads were subsampled with `seqtk sample` (well under 1% of
   the reads were retained).
2. The subsampled reads were classified with **Kraken2 version 2.17.1**
   against the **RefSeq complete genomes** database
   (`Kraken2_Bracken_RefSeq_Genomes_Complete`, downloaded in October 2024),
   in six separate runs, one per `--confidence` value listed above, each
   writing a Kraken report (`--report`).
3. Each Kraken report was converted to MPA format with KrakenTools
   `kreport2mpa.py`.
4. To keep the package small, the reports were reduced to a subset of the
   classified taxa. In the resulting files every row keeps its complete
   parent lineage, and taxa with a single read are present, i.e. the
   reduction was not a minimum read-count filter.

The exact commands (seqtk seed and fraction, Kraken2 options other than
`--confidence`, and the rule used to select taxa in step 4) were not
recorded. The files are provided to demonstrate and test the package, not as
a reproducible analysis. The general procedure for producing input for
karioCaS from new data is:

```bash
DB=/path/to/Kraken2_database
for CS in 0 0.2 0.4 0.6 0.8 0.9; do
    CODE=$(printf "%02d" "$(echo "$CS * 10" | bc | cut -d. -f1)")
    kraken2 --db "$DB" --confidence "$CS" --threads 8 \
        --report SAMPLE01_CS${CODE}.kreport --output /dev/null \
        reads_R1.fastq.gz reads_R2.fastq.gz   # add --paired for paired-end reads
    kreport2mpa.py -r SAMPLE01_CS${CODE}.kreport -o SAMPLE01_CS${CODE}.mpa
done
```

The resulting `SAMPLE01_CSxx.mpa` files go into `<project>/000_mpa_original/`
(see the package vignette).
