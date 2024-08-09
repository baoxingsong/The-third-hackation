import pandas as pd
import sys


df_node = pd.read_csv(sys.argv[1],sep=' ')
df_node.columns = ['gene', 'node']
with open('merge.order.txt', mode='w') as f, open('merge.order.chr.txt', mode='w') as f2:
    for i, arg in enumerate(sys.argv):
        if i < 2:
            continue
        df = pd.read_csv(sys.argv[i], sep='\t')
        df.columns = ['chr', 'gene', 'start', 'end']
        df = df.sort_values(by=['chr', 'start'], ascending=[True, True])
        df_merge = pd.merge(df, df_node, how='inner', on='gene')
        df_merge.to_csv(f'{sys.argv[i]}.order', index=False, sep='\t')
       
        name = f'{sys.argv[i]}'.split('.')[0]
        with open(f'{name}.all.sequence', mode='w') as f3:
            for chr,data in df_merge.groupby('chr'):
                node_list = [str(i) for i in data['node'].tolist()]
                node_str = ' '.join(node_list)
                f3.write(node_str)
                f3.write('\n')
        
        with open(f'{name}.all.sequence.genename', mode='w') as f4:
            for chr,data in df_merge.groupby('chr'):
                gene_list = [str(i) for i in data['gene'].tolist()]
                gene_str = ' '.join(gene_list)
                f4.write(gene_str)
                f4.write('\n')


        bed_name = f'{sys.argv[i]}'.split('.')[0]
        for chr,data in df_merge.groupby('chr'):
            node_list = [str(i) for i in data['node'].tolist()]
            f2.write(f'{bed_name}\t{chr}\n')
            node_str = ' '.join(node_list)
            f.write(node_str)
            f.write('\n')
