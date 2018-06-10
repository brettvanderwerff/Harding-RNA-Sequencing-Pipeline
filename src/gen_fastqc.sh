#!/bin/bash
fastq_dir=$1
fastqc_dir=$2
output_dir=$3
cd $fastq_dir
for fastqfile in *.fastq.gz
    do
        $fastqc_dir $fastqfile \
        --outdir=$output_dir
    done
