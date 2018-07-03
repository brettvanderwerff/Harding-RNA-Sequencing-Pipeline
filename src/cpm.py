import pandas as pd
from statistics import mean

def gen_cpm(count_matrix):
    '''Takes a gene count matrix as input and calculated the "counts per million" of each feild.
    '''
    feature_count_matrix = pd.read_csv(count_matrix, index_col='Geneid')
    cpm_df = pd.DataFrame(index=feature_count_matrix.index)
    for column in list(feature_count_matrix.columns.values):
        column_total = feature_count_matrix[column].sum()
        column_list = [(field/column_total)*1000000 for field in feature_count_matrix[column]]
        cpm_df[column] = pd.Series(column_list).values
    return cpm_df

def filter_matrix(cpm_df, count_matrix):
    '''Keeps only the rows of a dataframe that contain at least 2 fields with a "counts per million" value that
    are greater than .5
    '''
    feature_count_matrix = pd.read_csv(count_matrix, index_col='Geneid')
    cpm_df_nan = cpm_df[(cpm_df > 0.5)]
    cpm_df_nan_drop = cpm_df_nan.dropna(thresh=2)
    return feature_count_matrix.merge(cpm_df_nan_drop, on=list(feature_count_matrix.columns.values), left_index=True, right_index=True)

def tmm(filtered_matrix):
    '''Normalizes the filtered matrix values via the trimmed mean m values method
    see: 10.1186/gb-2010-11-3-r25. The reference sample is chosen via the EdgeR packages's method.
    '''
    column_sums = {}
    for column in list(filtered_matrix.columns.values):
        column_sums[column] = filtered_matrix[column].sum()
    scaled_df = pd.DataFrame(index=filtered_matrix.index)
    for column in list(filtered_matrix.columns.values):
        column_total = filtered_matrix[column].sum()
        column_list = [field / column_total for field in filtered_matrix[column]]
        scaled_df[column] = pd.Series(column_list).values
    quantiles = scaled_df.quantile(q=.75)
    average_quantile = mean(scaled_df.quantile(q=.75).values)

if __name__ == "__main__":
    count_matrix = r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\featureCounts_matrix.csv'
    cpm_df = gen_cpm(count_matrix)
    filtered_matrix = filter_matrix(cpm_df, count_matrix)
    tmm(filtered_matrix)






