import os

source_folder = '/media/vande060/101A04EA101A04EA/MSP_RNA_Seq/'

CONFIGURATION = {
    'genome_annotation_dir' : source_folder + '/genome/annotations/genes.gtf',
    'fastq_dir' : source_folder + 'fastqs',
    'fastqc_dir' : os.path.join(os.getcwd(), os.path.basename('FastQC'), os.path.basename('fastqc')),
    'fastqc_output_dir' : source_folder + 'FastQC_output',
    'genome_index_dir' : source_folder + 'genome_index',
    'genome_sequence_dir' : source_folder + '/genome/sequence/genome.fa',
    'bam_dir' : source_folder + 'bams',
    'star_dir' : os.path.join(os.getcwd(), os.path.basename('STAR-2.6.0a'), os.path.basename('bin'),
                            os.path.basename('Linux_x86_64'), os.path.basename('STAR'))
}