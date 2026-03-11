#!/bin/bash
set -eux

echo "Enter query FASTA file:"
read query

echo "Enter CARD BLAST database:"
read db

echo "Enter aro_index.tsv file:"
read aroindex

# Extract filename without extension
prefix=$(basename "$query")
prefix=${prefix%.*}

echo "Running BLAST search..."

blastn -query "$query" \
-db "$db" \
-evalue 1e-30 \
-out "${prefix}_blast_raw.tsv" \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

echo "Filtering hits with identity >= 90%..."

awk -F'\t' '$3 >= 90' "${prefix}_blast_raw.tsv" > "${prefix}_blast_filtered.tsv"

# Remove raw BLAST file to prevent accumulation
rm "${prefix}_blast_raw.tsv"

echo "Annotating ARO IDs..."

# Create output file with headers
echo -e "ARO_ID\tGene_Name\tAMR_Gene_Family\tDrug_Class\tResistance_Mechanism" > "${prefix}_resistance.tsv"

awk -F'\t' '
NR==FNR && FNR>1 {
    gene[$1]=$9
    drug[$1]=$10
    mech[$1]=$11
    next
}
{
    split($2,a,"|")
    aro=a[5]
    gene_name=a[6]

    if(aro in gene){
        print aro"\t"gene_name"\t"gene[aro]"\t"drug[aro]"\t"mech[aro]
    }
}
' "$aroindex" "${prefix}_blast_filtered.tsv" >> "${prefix}_resistance.tsv"

echo "Pipeline complete."

echo "Generated files:"
echo "${prefix}_blast_filtered.tsv"
echo "${prefix}_resistance.tsv"
