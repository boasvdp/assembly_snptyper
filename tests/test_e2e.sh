#!/bin/bash

set -euxo pipefail

TEST_DATA_DIR=tests/tmp

# Download 
mkdir -p "$TEST_DATA_DIR"
## M1UK strain from Lynskey et al.
curl --output "$TEST_DATA_DIR"/GCA_034818825.1.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_034818825.1/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
## M1global strain from Lynskey et al.
curl --output "$TEST_DATA_DIR"/GCA_034822885.1.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_034822885.1/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"

# Unzip
unzip -n "$TEST_DATA_DIR"/GCA_034818825.1.zip -d "$TEST_DATA_DIR"
unzip -n "$TEST_DATA_DIR"/GCA_034822885.1.zip -d "$TEST_DATA_DIR"

# List genomes
find "$TEST_DATA_DIR" -iname "*.fna" -print > "$TEST_DATA_DIR"/test_genomes.txt

# Install the script
python -m pip install .

# Run the script
assembly_snptyper \
    --list_input "$TEST_DATA_DIR"/test_genomes.txt \
    --vcf data/M1UK.vcf \
    --reference data/MGAS5005.fa \
    -p 2 > "$TEST_DATA_DIR"/test_output.txt

# Check the output
cmp "$TEST_DATA_DIR"/test_output.txt tests/test_expected_output.txt

rm -rf "$TEST_DATA_DIR"