import pandas as pd

gene_count_matrix = pd.read_csv(r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\trimmed_results.csv')

key = pd.read_csv(r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\key.csv')

example = pd.read_csv(r'C:\Users\vande060\Desktop\coding\projects\Harding-RNA-Sequencing-Pipeline\paper_result.csv')

genename = pd.read_csv('geneid_genename.txt', sep=',')

geneid_dd = genename.drop_duplicates(subset='Gene stable ID')

geneid_dd_renamed = geneid_dd.rename(columns={('Gene stable ID'): 'Geneid'})


merged = pd.merge(gene_count_matrix, geneid_dd_renamed[['Geneid', 'Gene name']], on='Geneid', how='left')

merge_key = pd.merge(merged, key, on='Geneid', how='left')
final = merge_key[merge_key.Sig != 0]

merge_final = pd.merge(final, example, on='Geneid', how='left')

merge_final.to_csv('final.csv', sep=',')

