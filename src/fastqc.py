import os
from src.config import CONFIGURATION
import subprocess
import wget
import zipfile


def get_fastqc(fastqc_dir):
    '''If FastQC is not already downloaded, function downloads FastQC version 0.11.7 zip file, unzips file, and
    makes fastqc executable.
    '''
    file_dir = os.path.dirname(__file__)
    if not os.path.isdir(os.path.join(file_dir, os.path.basename('FastQC'))):
        fastqc_url = 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip'
        wget.download(fastqc_url, out=file_dir)
        with zipfile.ZipFile(os.path.join(file_dir, os.path.basename('fastqc_v0.11.7.zip')), "r") as zip_obj:
            zip_obj.extractall(file_dir)
        os.remove(os.path.join(file_dir, os.path.basename('fastqc_v0.11.7.zip')))
        os.chmod(fastqc_dir, 0o755)

def exec_fastqc(fastq_dir, fastqc_dir, fastqc_output_dir):
    '''Runs FastQC to analyze RNA read quality
    '''
    file_dir = os.path.dirname(__file__)
    if not os.path.isdir(fastqc_output_dir):
        os.mkdir(fastqc_output_dir)
    os.chmod('{}/gen_fastqc.sh'.format(file_dir), 0o755)
    subprocess.call(['{}/gen_fastqc.sh'.format(file_dir), fastq_dir, fastqc_dir, fastqc_output_dir])

def run():
    '''
    Runs the fastqc script.
    '''
    fastqc_dir = CONFIGURATION['fastqc_dir']
    fastqc_output_dir = CONFIGURATION['fastqc_output_dir']
    fastq_dir = CONFIGURATION['fastq_dir']
    get_fastqc(fastqc_dir)
    exec_fastqc(fastq_dir, fastqc_dir, fastqc_output_dir)


if __name__ == '__main__':
    run()