from src.config import CONFIGURATION
import pandas as pd
import os


def get_gene_counts(count_output_dir):
    '''Reads each raw output csv from subread.py (one csv for each sample), outputs a single csv containing the raw
    counts for each sample and a column for gene id
    '''
    for file in os.listdir(count_output_dir):
        if file.endswith('txt'):
            df = pd.read_csv(os.path.join(count_output_dir, os.path.basename(file)), sep='\t', skiprows=1)
            geneid = df[['Geneid']]
            break
    for file in os.listdir(count_output_dir):
            if file.endswith('txt'):
                df = pd.read_csv(os.path.join(count_output_dir, os.path.basename(file)), sep='\t', skiprows=1)
                trimmed_df = df.drop(['Chr', 'Start', 'Length', 'End', 'Strand'], axis=1)
                geneid[trimmed_df.columns.values[1]] = trimmed_df[trimmed_df.columns.values[1]]
    geneid.to_csv(os.path.join(count_output_dir,
                               os.path.basename('featureCounts_matrix.csv')), sep=',', index=False)


if __name__ == '__main__':
    count_output_dir = CONFIGURATION['count_output_dir']
    get_gene_counts(count_output_dir)


