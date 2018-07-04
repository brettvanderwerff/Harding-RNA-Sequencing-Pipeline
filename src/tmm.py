import pandas as pd
from src import cpm
from statistics import mean
pd.set_option('precision',10)


def id_ref(filtered_matrix):
    '''Identifies the reference sample by which all other samples will be normalized against by the trimmed mean m
    values method.
    '''
    scaled_df = pd.DataFrame(index=filtered_matrix.index)
    for column in list(filtered_matrix.columns.values):
        column_total = filtered_matrix[column].sum()
        column_list = [field / column_total for field in filtered_matrix[column]]
        scaled_df[column] = pd.Series(column_list).values
    quantiles = scaled_df.quantile(q=.75)
    average_quantile = mean(scaled_df.quantile(q=.75).values)
    return (abs(quantiles - average_quantile)).idxmin()


def tmm(filtered_matrix):
    '''Normalizes the filtered matrix values via the trimmed mean m values method
    see: 10.1186/gb-2010-11-3-r25. The reference sample is chosen via the EdgeR packages's method.
    '''
    return None


if __name__ == '__main__':
    count_matrix = r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\featureCounts_matrix.csv'
    cpm_df = cpm.gen_cpm(count_matrix)
    filtered_matrix = cpm.filter_matrix(cpm_df, count_matrix)
    print(id_ref(filtered_matrix))