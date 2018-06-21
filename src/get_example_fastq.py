import os
import wget

def get_example_fastq():
    '''Downloads a set of example .fastqs for demonstration, only to be used if the user does not have their own .fastq
    files.
    '''
    out_dir = os.path.join(os.path.dirname(__file__),
                           os.path.basename('Example'),
                           os.path.basename('fastqs'))
    sample = {
        '001': 'SRR1036961',
        '002' : 'SRR1036962',
        '005' :'SRR1036965',
        '006' : 'SRR1036966',
        '007' : 'SRR1036967',
        '003' : 'SRR1036963'}
    for key, value in sample.items():
        wget.download('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/{}/{}/{}.fastq.gz'.format(key, value, value), out=out_dir)

def run():
    '''Runs the get_example_fastq script.
    '''
    get_example_fastq()

if __name__ == "__main__":
    get_example_fastq()

