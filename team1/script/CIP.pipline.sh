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


for i in "${array[@]}";do
	for j in "${array[@]}";do
		if [ "$i" = "$j" ];then continue;fi
   		$software_pwd/diamond blastp -d $j -q $data_pwd/$i.pep.fa -o "$i"_"$j".blastout -p $CPU -k 1 -e 1e-5 -f 6 --outfmt 6
	done
done

export PERL5LIB=/media/songlab/14t1/Team1/software/miniconda3/pkgs/perl-bioperl-core-1.007002-pl5262hdfd78af_3/lib/site_perl/5.26.2/:$PERL5LIB
for i in "${array[@]}";do
        for j in "${array[@]}";do
                if [ "$i" = "$j" ];then continue;fi
		cp "$i"_"$j".blastout "$i"_"$j".blastout.CIP_CALP-out
		# perl $script_pwd/6.CIP_CALP.pl "$i"_"$j".blastout > "$i"_"$j".blastout.CIP_CALP-out
	done
done

cd ..


# MCScanx
mkdir 02_MCScanx
cd 02_MCScanx
for i in "${array[@]}";do
        for j in "${array[@]}";do
                if [ "$i" = "$j" ];then
                $software_pwd/diamond blastp -d ../01_diamond/$i -q $data_pwd/$i.pep.fa -o "$i".blast -p $CPU -k 6 -e 1e-5 -f 6 --outfmt 6
        	cat $data_pwd/"$i".bed > "$i".gff
                MCScanX "$i"
		fi
	done
done
cd ..

# drop_tandem
mkdir 03_drop_tandem
cd 03_drop_tandem
cat ../01_diamond/*.CIP_CALP-out |grep -v '#' | sort |uniq > all.blast
cat ../02_MCScanx/*.tandem |sed 's/,.*//g' > all.tandem 
grep -v -F -f all.tandem all.blast > all.drop_tandem.blast
cut -f1,2 all.drop_tandem.blast |sed 1d |grep -v 'Query'> all.drop_tandem.blast.pair
python $script_pwd/cluster_nodes.py all.drop_tandem.blast.pair all.drop_tandem.blast.pair.cluster
for i in "${array[@]}";do
	cp $data_pwd/"$i".bed .
done
python $script_pwd/group_matrix.py
python $script_pwd/get_chr_order.py all.drop_tandem.blast.pair.cluster.filter `ls *.bed`
cd ..

# drimm.sequence
mkdir 04_DrimmSequence
cd 04_DrimmSequence
cp ../03_drop_tandem/merge.order.txt drimm.sequence
mono $script_pwd/Program.exe ./drimm.sequence ./ 20 $Drimmsynteny_WGDnum


mkdir processDrimm_result
cd processDrimm_result
cp ../../03_drop_tandem/*.sequence ../../03_drop_tandem/*.genename .
cd ..
cp $script_pwd/processDrimm.py .
cp -r $script_pwd/utils/ .
# sed -i "s/chr_number = .*/chr_number = $chr_number/g" processDrimm.py
# sed -i "s/sp_list = .*/sp_list = $sp_list/g" processDrimm.py
# sed -i "s/target_rate = .*/target_rate = \'$target_rate\'/g" processDrimm.py
python processDrimm.py ../03_drop_tandem/merge.order.chr.txt
cd processDrimm_result/finalBlocks
for i in `ls *.block|sed 's/\..*//g'`;do echo "> $i" >> mgra.block;sed 's/$/\$/g;s/^s //g' $i.final.block >> mgra.block;echo "" >> mgra.block;done
cd ..
cp $script_pwd/processGenenumber.py .
python processGenenumber.py ../../03_drop_tandem/merge.order.chr.txt
cd ../../


#mgra
mkdir 05_mgra
cd 05_mgra
cp ../04_DrimmSequence/processDrimm_result/finalBlocks/mgra.block  .

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
mkdir 06_plot
cd 06_plot
cp $script_pwd/chromosomeRearrangementPainting.py .
cp $script_pwd/BlockMatchingOptimization.py .
cp ../05_mgra/out/genomes/*.block .

an_file=$(find . -type f -name "*.gen.block" | awk -F/ '{print $NF, length($NF)}' | sort -k2 -nr | head -n 1 | awk '{print $1}')
for i in `ls *.gen.block`;do if [ "$i" = "$an_file" ]; then continue;else python chromosomeRearrangementPainting.py ../04_DrimmSequence/processDrimm_result/finalBlocks/blockindex.genenumber $i $an_file;fi; done
