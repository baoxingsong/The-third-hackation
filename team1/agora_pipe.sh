#!/usr/bin/bash
### all the step in conda envs agora conda activate agora
ortho="orthofinder/OrthoFinder/Results_Aug06/Phylogenetic_Hierarchical_Orthogroups/"
pre_dir="00.pre"
gene_dir="01.gene.list"
result_dir="02.result"

if [ ! -d $pre_dir ]; then mkdir $pre_dir; fi
if [ ! -d $gene_dir ]; then mkdir $gene_dir; fi
if [ ! -d $result_dir ]; then mkdir $result_dir; fi

### prepare the input of agora
python3 ../software/Agora/src/import/convert_hogs_sp.py -of_hogs ${ortho} -outdir ${pre_dir}
cp ${ortho}../Species_Tree/SpeciesTree_rooted_node_labels.txt ./Species.tree
ls sp/*.gff |sed 's/sp\///' |sed 's/.gff//' > sp_name

for sp in $(cat sp_name);
do
	awk -F '\t|;' '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' |sed "s/ID=/${sp}_/"|sed 's/+/1/' |sed 's/-/-1/' > ${gene_dir}/genes.${sp}.list
done

### running agora

../software/Agora/src/agora-plants.py ./Species.tree ${pre_dir}/orthologyGroups.%s.list ${gene_dir}/genes.%s.list -workingDir=${result_dir} -nbThreads=8

### visualize in 02.result

../software/Agora/src/misc.compareGenomes.py  ancGenomes/plants-workflow/ancGenome.N3.list.bz2 ancGenomes/plants-workflow/ancGenome.N1.list.bz2 ancGenes/size-0.9-1.1/ancGenes.N1.list.bz2 -mode=drawKaryotype -minChrSize=100 +sortBySize > N3.ps
