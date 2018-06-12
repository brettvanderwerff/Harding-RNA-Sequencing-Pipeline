import os
from src.config import CONFIGURATION
import subprocess
import wget
import zipfile

def get_star():
    '''
     If STAR is not already downloaded, function downloads STAR version 2.6.0a zip file, and unzips file.
     '''
    if not os.path.isdir('STAR-2.6.0a'):
        star_url = 'https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz'
        wget.download(star_url)
        with zipfile.ZipFile("2.6.0a.tar.gz", "r") as zip_obj:
            zip_obj.extractall()
        os.remove(os.path.join(os.getcwd(), os.path.basename('2.6.0a.tar.gz')))

def make_genome_index(genome_annotation_dir, genome_sequence_dir, genome_index_dir, star_dir):
    '''
    Checks if a genome index already exists in the 'genome index directory'. If none exists STAR is used to create one
    from the sequence and annotations of a genome.
    '''
    if not os.path.isfile(os.path.join(genome_index_dir, os.path.basename('genome'))):
        os.chmod('./gen_genome.sh', 0o755)
        subprocess.call(['./gen_genome.sh', genome_annotation_dir, genome_sequence_dir, genome_index_dir, star_dir])


def exec_map(genome_index_dir, fastq_dir, star_dir, bam_output_dir):
    '''Uses STAR to map RNA reads to a genome index
    '''
    os.chmod('./map_reads.sh', 0o755)
    subprocess.call(['./map_reads.sh', genome_index_dir, fastq_dir, star_dir, bam_output_dir])


if __name__ == '__main__':
    genome_annotation_dir = CONFIGURATION['genome_annotation_dir']
    genome_index_dir = CONFIGURATION['genome_index_dir']
    genome_sequence_dir = CONFIGURATION['genome_sequence_dir']
    fastq_dir = CONFIGURATION['fastq_dir']
    star_dir = CONFIGURATION['star_dir']
    bam_dir = CONFIGURATION['bam_dir']
    exec_map(genome_sequence_dir, fastq_dir, star_dir, bam_dir)