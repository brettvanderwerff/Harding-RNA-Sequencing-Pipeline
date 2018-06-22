# Harding-RNA-Sequencing-Pipeline
A pipeline for analyzing RNA sequencing data

==WIP==

## Purpose

The purpose of this tool is to provide a partially automated analysis of raw RNA sequencing data from a two condition experiment to produce a list of differentially expressed genes. The pipeline for this tool is as follows:

Single end .fastq file --> FASTQC --> STAR --> Subread --> EdgeR

The EdgeR analysis was inspired by the following publication: doi: 10.12688/f1000research.8987.2

## Dependencies 

*General Dependencies*

* Linux (tested on Ubuntu 18.04)
* Python 3 (tested on 3.6)
* R (tested on 3.4.3)
* 32 GB of RAM (STAR requires 32 GB of RAM for alignment to human genome) 
* Over 100 GB of free HD space 

*Python Dependencies*

* pandas==0.23.0
* wget==3.2

*R Dependencies*

* edgeR==3.20.9

## Setup and Operation

1. Clone repo to your computer
2. Open `src/config.py` and change the source_folder variable to point to the directory containing your project. Your project should have the following structure:

```
.
├── fastqs  # folder for raw .fastq files
├──genome
│   ├── sequence # folder to place a FASTA format genome sequence file (must be named 'genome.fa')
│   └── annotation # folder to place genome annotation (must be named 'genes.gtf')
└── genome_index  # folder STAR will write the genome index to
```

Note that by default the source_folder is set to an 'Example' folder to where .fastq files from the following study will be downloaded to: https://www.ebi.ac.uk/ena/data/view/PRJNA229803 for an example analysis. If using your own project, prevent the downloading of these files by commenting out the `get_example_fastq` function in the `rna_seq_analysis.py` module

Also note that a copy of the human genome in FASTA format (for placing in the 'sequence' folder) can be obtained here:

http://ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

And the corresponding gene annotations for the human genome (for placing in the 'annotation' folder) can be obtained here: 

http://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz

3. Open `src/diff_gene_expression.R` and change any reference to grouping factors to fit your study (Note that the grouping factors are 'msp_ron' and 'control' by default to match the example .fastq files)

4. Begin the analysis by running the `rna_seq_analysis.py` file, analysis may take anywhere from a few hours to a few days depending on the size of the genome and the number of .fastq files being analyzed. Results will be written to the 'featureCounts_output' directory of your project folder. 

           



