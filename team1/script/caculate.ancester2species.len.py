import pandas as pd
import sys

block_file = sys.argv[1]
genenumber_genename = sys.argv[2]
bed_file = sys.argv[3]

# python caculate.ancester2species.len.py ../../../06_plot/BCAE.gen.block Orin.final.synteny.genename /media/songlab/14t1/Team1/CIP/data/Orin.bed

block_list = []
with open(block_file, mode = 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split(' ')[1:]

        for i in line:
            if i.startswith('-'):
                i = i.replace('-','')
            block_list.append(str(i))

df_genename = pd.DataFrame(columns=['block', 'chr', 'direction', 'start_gene', 'end_gene'])
with open(genenumber_genename, mode = 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split(':')
        if line[3].split(' ')[0] == '-':
            df_genename.loc[len(df_genename.index)] = [line[0], line[2], '-', line[3].split(' ')[-1], line[3].split(' ')[1]]
        else:
            df_genename.loc[len(df_genename.index)] = [line[0], line[2], '+', line[3].split(' ')[1], line[3].split(' ')[-1]]
df_genename = df_genename[df_genename['block'].isin(block_list)]
df_bed =  pd.read_csv(bed_file, header=None, sep='\t')
df_bed.columns=['chr', 'gene', 'start_num', 'end_num']
df_genename = pd.merge(df_genename, df_bed[['gene', 'start_num', 'end_num']], how='left', left_on = 'start_gene', right_on = 'gene')
df_genename = pd.merge(df_genename, df_bed[['gene', 'start_num', 'end_num']], how='left', left_on = 'end_gene', right_on = 'gene', suffixes=('_start','_end'))
df_genename['block_length'] = df_genename['end_num_end'] - df_genename['start_num_start']
ancester_length = df_genename['block_length'].sum()
genome_length = 0
for i,data in df_bed.groupby('chr'):
    genome_length += data['end_num'].max()
print(ancester_length/genome_length)


