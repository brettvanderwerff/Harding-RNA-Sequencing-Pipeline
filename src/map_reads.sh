#!/bin/bash
genome_index_dir=$1
fastq_dir=$2
star_dir=$3
bam_output_dir=$4

cd $fastq_dir
for fastqfile in *.fastq.gz
do
    base=$(basename "$fastqfile" .fastq.gz)
    echo "mapping $base to genome"
	$star_dir \
		--genomeDir $genome_index_dir \
		--runThreadN 1 \
		--readFilesIn $fastqfile \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix $bam_output_dir$base
done
