#!/bin/bash

data_pwd='/media/songlab/14t1/Team1/CIP/data'
software_pwd='/media/songlab/14t1/Team1/software'
script_pwd='/media/songlab/14t1/Team1/CIP/script'
CPU=20
Drimmsynteny_WGDnum='11'
Tree='(((Seit,Sobi),(Orin,Hovu)),Anco)'
array=(Orin Seit Sobi Anco Hovu)

# DIAMOND
mkdir 01_diamond
cd 01_diamond

for i in "${array[@]}";do
   $software_pwd/diamond makedb -d $i --in $data_pwd/$i.pep.fa
done

# one2one
for i in "${array[@]}";do
	for j in "${array[@]}";do
		if [ "$i" = "$j" ];then continue;fi
		echo "${i}_${j}"
   		$software_pwd/diamond blastp -d $j -q $data_pwd/$i.pep.fa -o "$i"_"$j".blastout -p $CPU -k 6 -e 1e-5 -f 6 --outfmt 6
		python3 $script_pwd/prepareInput.py -q $data_pwd/$i.gff -r $data_pwd/$j.gff -b "$i"_"$j".blastout -o "$i"_"$j".prepare -qp "${i}_" -rp "${j}_"
		$software_pwd/AnchorWave/anchorwave pro -i "$i"_"$j".prepare -o "$i"_"$j".prepare.anchorwave -R 1 -Q 1
	done
done

cd ..


# drop_tandem
mkdir 02_drop_tandem
cd 02_drop_tandem
cat ../01_diamond/*.anchorwave |grep -v '#' | grep -v 'refGene'|cut -f1,6 | sort |uniq > all.anchorwave.pair
python $script_pwd/cluster_nodes.py all.anchorwave.pair all.anchorwave.pair.cluster
gff_wd=''
for i in "${array[@]}";do
    gff_wd+="${i}:${data_pwd}/${i}.gff;"
done
gff_wd="${gff_wd%?}"
# echo "python $script_pwd/cleanCluster.py -c all.anchorwave.pair.cluster -o all.anchorwave.pair.cluster.filter -g "${gff_wd}" -m 5"
python $script_pwd/cleanCluster.py -c all.anchorwave.pair.cluster -o all.anchorwave.pair.cluster.filter -g "${gff_wd}" -m 5
for i in "${array[@]}";do
	cp $data_pwd/"$i".bed .
done
sed -i 's/\t/ /g' all.anchorwave.pair.cluster.filter
python $script_pwd/get_chr_order.py all.anchorwave.pair.cluster.filter `ls *.bed`
cd ..

# drimm.sequence
mkdir 03_DrimmSequence
cd 03_DrimmSequence
cp ../02_drop_tandem/merge.order.txt drimm.sequence
mono $script_pwd/Program.exe ./drimm.sequence ./ 20 $Drimmsynteny_WGDnum


mkdir processDrimm_result
cd processDrimm_result
cp ../../02_drop_tandem/*.sequence ../../02_drop_tandem/*.genename .
cd ..
cp $script_pwd/processDrimm.py .
cp -r $script_pwd/utils/ .
# sed -i "s/chr_number = .*/chr_number = $chr_number/g" processDrimm.py
# sed -i "s/sp_list = .*/sp_list = $sp_list/g" processDrimm.py
# sed -i "s/target_rate = .*/target_rate = \'$target_rate\'/g" processDrimm.py
python processDrimm.py ../02_drop_tandem/merge.order.chr.txt
cd processDrimm_result/finalBlocks
for i in `ls *.block|sed 's/\..*//g'`;do echo "> $i" >> mgra.block;sed 's/$/\$/g;s/^s //g' $i.final.block >> mgra.block;echo "" >> mgra.block;done
cd ..
cp $script_pwd/processGenenumber.py .
python processGenenumber.py ../../02_drop_tandem/merge.order.chr.txt
cd ../../


#mgra
mkdir 04_mgra
cd 04_mgra
cp ../03_DrimmSequence/processDrimm_result/finalBlocks/mgra.block  .

English_chr=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)
echo '[Genomes]' > cfg.file
array_length=${#array[@]}
echo $Tree > Tree.file
for (( i=0; i<array_length; i++ ))
do
	echo "${English_chr[$i]} ${array[$i]}" >> cfg.file
	chr1=${English_chr[$i]}
	chr2=${array[$i]}
	sed -i "s/"$chr2"/"$chr1"/g" Tree.file
	sed -i "s/"$chr2"/"$chr1"/g" mgra.block
done
echo "" >> cfg.file
echo "[Trees]" >> cfg.file
cat Tree.file >> cfg.file
rm Tree.file

mkdir out
mgra -c cfg.file -g mgra.block -o out
cd out/genomes
for i in `ls *.gen|sed 's/\..*//g'`;do grep -v '#' $i.gen |grep -v '^$'|sed 's/^/s /g;s/\$//g;s/+//g' > $i.gen.block;done
cd ../../..

# plot
mkdir 05_plot
cd 05_plot
cp $script_pwd/chromosomeRearrangementPainting.py .
cp $script_pwd/BlockMatchingOptimization.py .
cp ../04_mgra/out/genomes/*.block .

an_file=$(find . -type f -name "*.gen.block" | awk -F/ '{print $NF, length($NF)}' | sort -k2 -nr | head -n 1 | awk '{print $1}')
for i in `ls *.gen.block`;do if [ "$i" = "$an_file" ]; then continue;else python chromosomeRearrangementPainting.py ../03_DrimmSequence/processDrimm_result/finalBlocks/blockindex.genenumber $i $an_file;fi; done


