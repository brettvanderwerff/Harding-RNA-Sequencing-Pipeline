import pandas as pd

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

def filter_cpm(cpm_df):
    cpm_df_nan = cpm_df[(cpm_df >= 0.5)]
    cpm_df_nan_drop = cpm_df_nan.dropna(thresh=2)
    return cpm_df.merge(cpm_df_nan_drop, on='Geneid')

if __name__ == "__main__":
    cpm_df = gen_cpm(count_matrix='featureCounts_matrix.csv')
    filtered_df = filter_cpm(cpm_df)
    print(filtered_df.head().to_string())


