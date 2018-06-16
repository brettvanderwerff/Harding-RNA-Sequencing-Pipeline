import os
import subprocess
from src.config import CONFIGURATION
import wget
import zipfile

def get_subread():
    '''If Subread is not already downloaded, function downloads Subread version 1.6.2 zip file and unzips file.
    '''
    if not os.path.isdir(os.path.join(os.path.dirname(__file__), os.path.basename('subread-1.6.2-Linux-x86_64'))):
        subread_url = 'https://www.sourceforge.net/projects/subread/files/subread-1.6.2-Linux-x86_64.tar.gz/download'
        wget.download(subread_url)
        with zipfile.ZipFile("subread-1.6.2-Linux-x86_64.tar.gz", "r") as zip_obj:
            zip_obj.extractall()


def exec_count(feature_counts_dir, bams_dir, annotations_dir, count_output_dir):
    '''Uses the featureCounts program of Subread to count the number of reads that align with specific genes within
    the genome.
    '''
    if not os.path.isdir(count_output_dir):
        os.mkdir(count_output_dir)
    os.chmod('./feature_counts.sh', 0o755)
    subprocess.call(['./feature_counts.sh', feature_counts_dir, bams_dir, annotations_dir, count_output_dir])

def run():
    '''Runs the subread script.
    '''
    bam_output_dir = CONFIGURATION['bam_output_dir']
    feature_counts_dir = CONFIGURATION['feature_counts_dir']
    genome_annotation_dir = CONFIGURATION['genome_annotation_dir']
    count_output_dir = CONFIGURATION['count_output_dir']
    get_subread()
    exec_count(feature_counts_dir, bam_output_dir, genome_annotation_dir, count_output_dir)


if __name__ == '__main__':
    run()