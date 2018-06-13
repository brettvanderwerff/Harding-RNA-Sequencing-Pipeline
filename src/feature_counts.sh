#!/bin/bash
feature_counts_dir=$1
bams_dir=$2
annotations_dir=$3
output_dir=$4

cd $bams_dir
pwd
for bamfile in *Aligned.out.bam
do
    base=$(basename "$bamfile" Aligned.out.bam)
    echo "Counting features of $base"
	$feature_counts_dir \
	    -T 6 \
	    -t gene \
	    -g gene_id \
	    -a $annotations_dir \
	    -o $output_dir$base.count.txt \
	    -s 0 \
	    $bamfile

done
