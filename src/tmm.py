from math import log
import numpy as np
import pandas as pd
from src import cpm
from statistics import mean
pd.set_option('precision',10)

def scale_df(filtered_matrix):
    '''Finds the scaled counts for each gene by dividing the read counts for each gene by the total number of counts
    in the library
    '''
    scaled_df = pd.DataFrame(index=filtered_matrix.index)
    for column in list(filtered_matrix.columns.values):
        column_total = filtered_matrix[column].sum()
        column_list = [field / column_total for field in filtered_matrix[column]]
        scaled_df[column] = pd.Series(column_list).values
    return scaled_df

def id_ref(scaled_df):
    '''Identifies the reference sample by which all other samples will be normalized against by the trimmed mean m
    values method.
    '''
    quantiles = scaled_df.quantile(q=.75)
    average_quantile = mean(scaled_df.quantile(q=.75).values)
    return (abs(quantiles - average_quantile)).idxmin()

def log_ratio(reference_sample, scaled_df):
    '''Finds the ratio of each gene count divided to total library size for each sample. This ratio for each gene in
    each sample is then divided by the corresponding ratio of the refrence sample. The log2 of this division is then
    performed.
    See http://evowiki.haifa.ac.il/images/f/f5/A_scaling_normalization_method_for_differential_expression_%28TMM%29.pdf
    '''
    log_ratio_df = pd.DataFrame(index=scaled_df.index)
    for column in list(scaled_df.columns.values):
        column_list = []
        for field_subject, field_reference in zip(scaled_df[column], scaled_df[reference_sample]):
            try:
                column_list.append(log(field_subject,2)/log(field_reference,2))
            except ValueError as e:
                column_list.append(np.nan)
        log_ratio_df[column] = pd.Series(column_list).values
    return log_ratio_df.dropna()


def tmm(filtered_matrix):
    '''Normalizes the filtered matrix values via the trimmed mean m values method
    see: 10.1186/gb-2010-11-3-r25. The reference sample is chosen via the EdgeR packages's method.
    '''
    return None


if __name__ == '__main__':
    count_matrix = r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\featureCounts_matrix.csv'
    cpm_df = cpm.gen_cpm(count_matrix)
    filtered_matrix = cpm.filter_matrix(cpm_df, count_matrix)
    scaled_df = scale_df(filtered_matrix)
    reference_sample = id_ref(filtered_matrix)
    log_ratio(reference_sample, scaled_df)
