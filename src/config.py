import os

source_folder = '/media/vande060/101A04EA101A04EA/PDGF_FGF_4_hr/'

CONFIGURATION = {
    'count_output_dir' : source_folder + 'featureCounts_output/',
    'genome_annotation_dir' : source_folder + 'genome/annotation/genes.gtf',
    'fastq_dir' : source_folder + 'fastqs',
    'fastqc_dir' : os.path.join(os.path.dirname(__file__),
                                os.path.basename('FastQC'),
                                os.path.basename('fastqc')),
    'fastqc_output_dir' : source_folder + 'FastQC_output',
    'feature_counts_dir' : os.path.join(os.path.dirname(__file__),
                                        os.path.basename('subread-1.6.2-Linux-x86_64'),
                                        os.path.basename('bin'),
                                        os.path.basename('featureCounts')),
    'genome_index_dir' : source_folder + 'genome_index/',
    'genome_sequence_dir' : source_folder + 'genome/sequence/genome.fa',
    'bam_output_dir' : source_folder + 'bams/',
    'star_dir' : os.path.join(os.path.dirname(__file__),
                              os.path.basename('STAR-2.6.0a'),
                              os.path.basename('bin'),
                              os.path.basename('Linux_x86_64'),
                              os.path.basename('STAR'))
}