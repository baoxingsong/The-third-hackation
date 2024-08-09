import pandas as pd
import sys

# python caculate.ancester.num.py ../../../06_plot/BCAE.gen.block blockindex.genenumber 

block_file = sys.argv[1]
genenumber_file = sys.argv[2]
block_list = []
with open(block_file, mode = 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split(' ')[1:]
            
        for i in line:
            if i.startswith('-'):
                i = i.replace('-','')
            block_list.append(int(i))
df = pd.read_csv(genenumber_file, header = 0,sep='\t')
gene_num = df[df['blockID'].isin(block_list)]['blockLength'].sum()
print(gene_num)

