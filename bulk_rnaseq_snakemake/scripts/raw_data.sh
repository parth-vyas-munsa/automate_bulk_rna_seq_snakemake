#!/bin/bash

set -euo pipefail

PROJECT="/home/bioinformatics/bulk_rnaseq_snakemake"

SRA_DIR="$PROJECT/data/sra"
FASTQ_DIR="$PROJECT/data/fastq"
REFERENCE_DIR="$PROJECT/reference"
INDEX_DIR="$PROJECT/index"

ACCESSION_FILE="$PROJECT/accession_list.txt"

THREADS=4

mkdir -p "$SRA_DIR"
mkdir -p "$FASTQ_DIR"
mkdir -p "$REFERENCE_DIR"
mkdir -p "$INDEX_DIR"

if [ ! -f "$ACCESSION_FILE" ]; then
    echo "ERROR: accession_list.txt not found:"
    echo "$ACCESSION_FILE"
    exit 1
fi

# DOWNLOAD SRA + CONVERT TO FASTQ

while read -r SAMPLE
do

    # Skip empty lines
    if [ -z "$SAMPLE" ]; then
        continue
    fi

    echo
    echo "Processing: $SAMPLE"


 echo "Downloading $SAMPLE"

    prefetch "$SAMPLE" \
        --output-directory "$SRA_DIR"

 SRA_FILE="$SRA_DIR/$SAMPLE/$SAMPLE.sra"

    if [ ! -f "$SRA_FILE" ]; then
        echo "ERROR: SRA file not found:"
        echo "$SRA_FILE"
        exit 1
    fi

  # Convert SRA → FASTQ

echo "Converting $SAMPLE to FASTQ"

    fasterq-dump \
        "$SRA_FILE" \
        --split-files \
        --threads "$THREADS" \
        --outdir "$FASTQ_DIR"


done < "$ACCESSION_FILE"

echo "COMPRESSING FASTQ FILES"

for FILE in "$FASTQ_DIR"/*.fastq
do

    if [ -f "$FILE" ]; then

        echo "Compressing: $FILE"

        gzip "$FILE"

    fi

done

# SHOW FASTQ FILES

echo
echo "FASTQ FILES:"
ls -lh "$FASTQ_DIR"

# DOWNLOAD REFERENCE GENOME

echo "DOWNLOADING MOUSE GRCm39 REFERENCE"

GENOME_URL="https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

GTF_URL="https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"


GENOME_GZ="$REFERENCE_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

GTF_GZ="$REFERENCE_DIR/Mus_musculus.GRCm39.116.gtf.gz"

GENOME="$REFERENCE_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"

GTF="$REFERENCE_DIR/Mus_musculus.GRCm39.116.gtf"


if [ ! -f "$GENOME_GZ" ] && [ ! -f "$GENOME" ]; then

    echo "Downloading GRCm39 genome..."

    wget -c \
        -O "$GENOME_GZ" \
        "$GENOME_URL"

else

    echo "Genome already downloaded."

fi


# DOWNLOAD GTF


if [ ! -f "$GTF_GZ" ] && [ ! -f "$GTF" ]; then

    echo "Downloading GTF annotation..."

    wget -c \
        -O "$GTF_GZ" \
        "$GTF_URL"

else

    echo "GTF already downloaded."

fi

if [ ! -f "$GENOME" ]; then

    echo "Decompressing genome..."

    gunzip "$GENOME_GZ"

fi

if [ ! -f "$GTF" ]; then

    echo "Decompressing GTF..."

    gunzip "$GTF_GZ"

fi

echo "BUILDING HISAT2 INDEX"

if [ ! -f "$INDEX_DIR/grcm39.1.ht2" ]; then

    echo "Building HISAT2 index..."

    hisat2-build \
        -p "$THREADS" \
        "$GENOME" \
        "$INDEX_DIR/grcm39"

else

    echo "HISAT2 index already exists."
    echo "Skipping index construction."

fi

echo "RAW DATA + REFERENCE SETUP COMPLETE"


echo "Project:"
echo "$PROJECT"


echo
echo "SRA directory:"
echo "$SRA_DIR"


echo
echo "FASTQ directory:"
echo "$FASTQ_DIR"


echo
echo "Reference directory:"
echo "$REFERENCE_DIR"


echo
echo "HISAT2 index directory:"
echo "$INDEX_DIR"


echo
echo "HISAT2 index files:"
ls -lh "$INDEX_DIR"


echo
echo "DONE"