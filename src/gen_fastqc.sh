#!/bin/bash
fastq_dir=$1
fastqc_dir=$2
fastqc_output_dir=$3
cd $fastq_dir
for fastqfile in *.fastq.gz
    do
        $fastqc_dir $fastqfile \
        --outdir=$fastqc_output_dir
    done
