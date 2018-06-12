#!/bin/bash
genome_annotation_dir=$1
genome_sequence_dir=$2
genome_index_dir=$3
star_dir=$4

$star_dir \
    --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir $genome_index_dir \
    --genomeFastaFiles $genome_sequence_dir \
    --sjdbGTFfile $genome_annotation_dir \
    --sjdbOverhang 100 \
    --genomeSAsparseD 2 \
    --limitIObufferSize 80000000

