#!/usr/bin/bash
### software install if already install it ignore
### conda create -n iags_auto python=3.9
### conda activate iags_auto
### conda install -c gurobi -c conda-forge -c huntguo iags_auto
### sudo apt install mono-devel
### Please attention that iags_auto need a license file gurobi.lic, the license file can get by following this markdown
### https://github.com/xjtu-omics/IAGS_AUTO/blob/main/utils/static/gurobi.md

### Chromosome ID in gff file must start with chr and the number should be 1-9
ortho=$(which orthofinder)
pep_fdr="01.pep"
ortho_result="01.pep/Orthofinder/Result_*/"
iags_dir="02.iags_auto"

if [ ! -d $pep_fdr ]; then mkdir $pep_fdr; fi
if [ ! -d $iags_dir ]; then mkdir $iags_dir; fi
### prepare the gff file and Orthogroup gene file
ls sp/*.gff3 |sed 's/sp\///'|sed 's/\.gff3//' > sp_name
for sp i in $(cat sp_name);
do
	awk -F '\t|;' '$3=="mRNA"{print "'${sp}'_"$1"\t"$9"\t"$4"\t"$5}' sp/${sp}.gff3 |sed "s/ID=/${sp}/" > ${iags_dir}/${sp}.gff
done

$ortho -f 01.pep/ -S diamond -M msa -T fasttree -t 24

cp ${ortho_result}/Orthogroups/Ortholog.tsv ${iags_dir}
sed "s/[0-9]//g" ${ortho_result}/Species_Tree/SpeciesTree_rooted.txt |sed "s/:\.//" > ${iags_dir}/species.tree
### use conda env to run iags_auto programm
#conda activate iags_auto
iags_auto -f ${iags_dir} -c 20 -d 13

#Result/painting is the directory of Ancestral code karyotype plot
