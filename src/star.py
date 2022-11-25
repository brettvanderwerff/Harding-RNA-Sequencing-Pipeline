import os
from src.config import CONFIGURATION
import subprocess
import wget
import tarfile

def get_star():
    '''
     If STAR is not already downloaded, function downloads STAR version 2.6.0a tars file, and extracts the tarfile.
     '''
    file_dir = os.path.dirname(__file__)
    if not os.path.isdir(os.path.join(file_dir, os.path.basename('STAR-2.6.0a'))):
        star_url = 'https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz'
        wget.download(star_url)
        with tarfile.open(os.path.join(file_dir, os.path.basename('"STAR-2.6.0a.tar.gz"')), "r") as tar_obj:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar_obj, file_dir)
        os.remove(os.path.join(file_dir, os.path.basename('STAR-2.6.0a.tar.gz')))

def make_genome_index(genome_annotation_dir, genome_sequence_dir, genome_index_dir, star_dir):
    '''
    Checks if a genome index already exists in the 'genome index directory'. If none exists STAR is used to create one
    from the sequence and annotations of a genome.
    '''
    file_dir = os.path.dirname(__file__)
    if not os.path.isdir(genome_index_dir):
        os.mkdir(genome_index_dir)
    if not os.path.isfile(os.path.join(genome_index_dir, os.path.basename('genome'))):
        os.chmod('{}/gen_genome.sh'.format(file_dir), 0o755)
        subprocess.call(['{}/gen_genome.sh'.format(file_dir), genome_annotation_dir, genome_sequence_dir, genome_index_dir, star_dir])


def exec_map(genome_index_dir, fastq_dir, star_dir, bam_output_dir):
    '''Uses STAR to map RNA reads to a genome index.
    '''
    file_dir = os.path.dirname(__file__)
    if not os.path.isdir(bam_output_dir):
        os.mkdir(bam_output_dir)
    os.chmod('{}/map_reads.sh'.format(file_dir), 0o755)
    subprocess.call(['{}/map_reads.sh'.format(file_dir), genome_index_dir, fastq_dir, star_dir, bam_output_dir])

def run():
    '''
    Runs the star script.
    '''
    genome_annotation_dir = CONFIGURATION['genome_annotation_dir']
    genome_index_dir = CONFIGURATION['genome_index_dir']
    genome_sequence_dir = CONFIGURATION['genome_sequence_dir']
    fastq_dir = CONFIGURATION['fastq_dir']
    star_dir = CONFIGURATION['star_dir']
    bam_output_dir = CONFIGURATION['bam_output_dir']
    get_star()
    make_genome_index(genome_annotation_dir, genome_sequence_dir, genome_index_dir, star_dir)
    exec_map(genome_index_dir, fastq_dir, star_dir, bam_output_dir)


if __name__ == '__main__':
    run()