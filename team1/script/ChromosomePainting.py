import sys

from util.chromosomeRearrangementPainting import plotChrsRearrangement


path = r'.'
path_fin = '../05_mgra/genomes'

# colors for target species chromosomes, preWGD yeast ancestor
colorlist = ['#DF1159','#1E93C9','#26AF67','#D5A1C5','#EBCA6D',
             '#94B51E','#000000','#E066FF']

block_length_file = path + '/processDrimm_result/finalBlocks/blockindex.genenumber'
target_species_block_file = path_fin + r'/EB.gen'
target_species_name = 'EB'
target_species_copy_number = 1

"""
preWGD yeast -> Fves
"""
rearranged_species_block_file = path_fin + r'/E.gen'
rearranged_species_name = 'E'
rearranged_species_copy_number = 1

outdir = path + '/plot'
plotChrsRearrangement(block_length_file,
                          rearranged_species_block_file,rearranged_species_name,rearranged_species_copy_number,
                          target_species_block_file,target_species_name,target_species_copy_number,
                          colorlist,outdir)
