import pandas as pd


df = pd.read_csv('all.drop_tandem.blast.pair.cluster', header=None, sep=' ')
df.columns = ['gene', 'group']

df['spe'] = df['gene'].str.split('_').str[0]
all_spe = set(df['spe'].to_list())
gene_matrix = pd.DataFrame(columns=list(all_spe)+ ['group'])
spe_dic = {}
for i,data in df.groupby('group'):
    spe_dic = {}
    for j in data.index:
        gene = data.loc[j,'gene']
        spe = data.loc[j,'spe']
        if spe not in spe_dic.keys():
            spe_dic[spe] = []
        spe_dic[spe].append(gene)
    
    matrix_num = len(gene_matrix.index)
    gene_matrix.loc[matrix_num,'group'] = i
    for j in list(all_spe):
        if j in spe_dic.keys():
            gene_matrix.loc[matrix_num,j] =  len(spe_dic[j])

gene_matrix = gene_matrix.fillna(0)
gene_matrix.to_csv('all.drop_tandem.blast.pair.cluster.genecount.csv', index=False)
select_df = gene_matrix[(gene_matrix[list(all_spe)] <= 1).all(axis=1)]
select_list = select_df['group'].to_list()
df_out = df[df['group'].isin(select_list)][['gene', 'group']]
df_out.to_csv('all.drop_tandem.blast.pair.cluster.filter', index=False,sep=' ')
