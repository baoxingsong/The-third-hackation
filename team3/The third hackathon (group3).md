# The third hackathon (group3)

## pangenome construction

### PGGB
Converting chrom to [PanSN-spec](https://github.com/pangenome/PanSN-spec) format

```
cat sample.list | while read f;
do 
	python3 simple_chr.py ${f}.fasta > ${f}.chr.fasta
	python3 fa-rename.py --ids change.list ${f}.chr.fasta > ${f}.newchr.fasta
	seqkit seq -u ${f}.newchr.fasta > ${f}.upper.fasta
done
```

```
>simple_chr.py
#!/usr/bin/python

import os
import sys

args=sys.argv
CDS={}
id=''
sequence=''

f=open(args[1],'r')
for line in f:
    if line.startswith('>C'):
        line = line.replace('\n', '')
        id = ''.join(line.replace('>','').split(' ')[6]).replace(',','')
        print('>Chr',id,sep = '')
    elif line.startswith('>J'):
        pass
    else:
        print(line,end = '')


f.close()
```



Runing PGGB

```
#!/bin/bash
#BATCH -o job.%j.out
#SBATCH -p tcum256c128Partition
#SBATCH -J job_name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
ulimit -u 40960
cat ${f}.upper.fasta > all.fasta
bgzip -@ 50 -c all.fasta > all.fasta.gz
samtools faidx all.fasta.gz
pggb -i all.fasta -o output2 -t 80 -T 80 -V GCA_036927085.1:# -M -m
```

### minigraph-cactus
```
#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p tcum256c128Partition
#SBATCH -J test_cactus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

source /data/sux/hackathon_group3/01.software/cactus-bin-v2.9.0/venv-cactus-v2.9.0/bin/activate
export PATH=$PATH:/home/sux/my_data/hackathon_group3/01.software/cactus-bin-v2.9.0/bin
export PATH=$PATH:/data/sux/hackathon_group3/01.software/cactus-bin-v2.9.0/venv-cactus-v2.9.0/bin

function copy(){  ##copy raw data to current directory ./data
for x in `ls ../00.data/*fna `;
do echo ${x} ;
cp ${x} ./data;
done
}
function RmcontigRenameHeader(){  ##Filter the contig/chloroplast seq in the copy data and rename the seq in the easiest character like '>Chr1'
for x in `ls ./data/*fna`;do
echo ${x}
seqkit seq ${x} -n|egrep 'contig|chloroplast'|cut -f1 > ${x}.rmseq;
seqkit grep -v -f ${x}.rmseq ${x} -n > ${x}.new.fa;
python rename.py ${x}.new.fa > ${x}.new.rename.fa
done
}
function rename(){  ##Rename the .fa files and generate the data.path file for the cactus input  
for x in `ls ./data/*new.rename.fa`;do
newname=`echo ${x}|cut -f 3 -d'/'|cut -f1 -d'.'`
mv ${x} ./data/${newname}.fa
cat ./data/${newname}.fa|grep '>'
done
}
function cactus(){
date
cactus-pangenome ./tmp ./data.path --outDir ./result --outName cactus_test --reference GCA_036927085 --consCores 20 --odgi --vcf --giraffe --gfa --gbz 
date
}
#cactus

```

## TE derived SV identification

>>**Stragety**
> ExampleData: Ri, Shahdara, Tsu, Yo, Oy(Ref) 

>1. Generate 4 samples' gvcf with the Oy reference genome based on **AnchorWave/Tassel**
>2. Generate 5 samples' TE annotation based on EDTA
>3. Convert vcf to inser.SV.bed and del.SV.bed seperately
>4. Convert TE to bed format 
>5. Use bedtools to intersect del.SV.bed with Oy.TEanno.bed
and inser.SV.bed with TE.anno.bed respectively.
>6. Extract the variants from the initial vcf files.

## What we do？
We write a python script named as 'TE_derived_SV.py' which generates the each samples' TE and SV bed files and embeds the **bedtools** software to get the overlap region and further obtain the variants coordinates.

```
python TE_derived_SV.py
cat ./result/*inser_TE_0.9.bed |awk '{print$6"\t"$11}'|sort -k1n,2n|uniq > ./result/all.inser.txt #Inser Alt genome
cat ./result/*del_TE_0.9.bed |awk '{print$6"\t"$7+1}'|sort -k1n,2n|uniq  > ./result/all.del.txt #Del Ref genome
```

```
>TEderivedSV.py
import sys,os
os.system('mkdir result')
def vcf_gff(vcf,gff,head):
	###Handle vcf file
	if vcf != "":
		print(head,vcf,'SV.bed building')
		wri_inser, wri_del = open(f"./result/{head}.inser.bed",'w'), open(f"./result/{head}.del.bed",'w')
		with open(vcf,'r') as f:
			for l in f:
				if '##' in l:continue
				if '#CHROM' in l:
					sample=l.strip().split()[9]
				else:
					l0=l.strip().split()
					ch,pos,ref, alt,queryinfo = l0[0],l0[1],l0[3],l0[4], l0[7].split(';')
					querych, querysta, queryEnd, querystrand = queryinfo[0].replace('ASM_Chr=',''), queryinfo[2].replace('ASM_Start=',''), queryinfo[1].replace('ASM_End=',''),queryinfo[3].replace('ASM_Strand=','')
					length = len(ref)-1
					querylength = abs(int(queryEnd)-int(querysta)+1)
					if querylength > 50: # Insertion
						if querystrand == "+":
							pri = f"{querych}\t{int(querysta)-1}\t{queryEnd}\tInser:{querylength}\t{ch}\t{pos}\t{int(pos)+length}\t{querystrand}\t{ref}\t{alt}"
						elif querystrand == "-":
							pri = f"{querych}\t{int(queryEnd)-1}\t{querysta}\tInser:{querylength}\t{ch}\t{pos}\t{int(pos)+length}\t{querystrand}\t{ref}\t{alt}"
						wri_inser.write(f'{pri}\n')
						#print(pri)
					if length+1 > 50: #Deletion
						pri = f"{ch}\t{int(pos)-1}\t{int(pos)+length}\tDel:{length+1}\t{querych}\t{int(querysta)}\t{queryEnd}\t{querystrand}\t{ref}\t{alt}"
						wri_del.write(f'{pri}\n')
						#print(pri)
		wri_inser.close()
		wri_del.close()
	###
	if gff != "":
		print(head,'TE.bed building')
		wri_te = open(f"./result/{head}.te.bed",'w')
		with open(gff,'r') as f:
			for l in f:
				if '##' in l:continue
				l0=l.strip().split()
				ch,sta,end,TEtype,strand=l0[0][3:],l0[3],l0[4],l0[2],l0[6]
				pri = f"{ch}\t{int(sta)-1}\t{end}\t{TEtype}\t{strand}"
				wri_te.write(f"{pri}\n")
	if vcf != "" and gff != "":
		print(f'TE derived SV processing')
		cmd = f"bedtools intersect -a ./result/{head}.te.bed -b ./result/{head}.inser.bed -F 0.9 -wa -wb \>./result/{head}.inser_TE_0.9.bed"
		os.system(f"echo {cmd} >> tmp.sh")
		cmd = f"bedtools intersect -a  ./result/Oy-0.te.bed  -b ./result/{head}.del.bed -F 0.9 -wa -wb \> ./result/{head}.del_TE_0.9.bed"
		os.system(f"echo {cmd} >> tmp.sh")	
		#os.system(f"bedtools intersect -a  ./result/Oy-0.te.bed  -b ./result/{head}.del.bed -F 0.9 -wa -wb> ./result/{head}.del_TE_0.9.bed")
		#os.system(f"bash tmp.sh")
		print("TE derived SV processing Done")
	print(f"\n")

wdf=os.listdir('.')
dic_wdf={}
vcf_gff('','Oy-0.mod.EDTA.TEanno.gff3','Oy-0')
for x in wdf:
	if 'gff' in x or 'vcf' in x:
		k=x.split('.')[0]
		if k not in list(dic_wdf):
			dic_wdf[k] = {}
		if 'gff' in x:
			dic_wdf[k]['gff'] = x
			
		elif 'vcf' in x:
			dic_wdf[k]['vcf'] = x
for x in list(dic_wdf):
	if 'gff' in list(dic_wdf[x]) and 'vcf' in list(dic_wdf[x]): 
		gff, vcf = dic_wdf[x]['gff'], dic_wdf[x]["vcf"]
		print('InputGff:',gff,'\tInputVcf:',vcf,'\tSample:',x)
		vcf_gff(vcf,gff,x)
```


## Inversion and translocation extraction

Running AnchorWave using Oy-0 as reference and Yo-0 as query

```
#!/bin/bash
#SBATCH -o job.anchorwave_.%j.out
#SBATCH -p tcum256c128Partition
#SBATCH -J TWGG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
echo Start Time is `date`
time1=`date +"%Y-%m-%d %H:%M:%S"`
#---------------------------------------------------------------------------
###ref is Oy-0
anchorwave gff2seq -r /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Oy-0.fasta -i /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Oy-0.gff -o Oy-0_cds.fa
minimap2 -x splice -t 15 -k 12 -a -p 0.4 -N 20 --cs /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Oy-0.fasta Oy-0_cds.fa > Oy-0.sam
minimap2 -x splice -t 50 -k 12 -a -p 0.4 -N 20 --cs /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Yo-0.fasta Oy-0_cds.fa > Yo-0.sam
anchorwave proali -r /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Oy-0.fasta -i /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Oy-0.gff -ar Oy-0.sam -as Oy-0_cds.fa -a Yo-0.sam -s /home/sux/my_data/hackathon_group3/anchorwave_vg/00.raw_data/Yo-0.fasta -n Oy-0_Yo-0.anchor -o Oy-0_Yo-0.maf -f Oy-0_Yo-0.fragment.maf -R 1 -Q 1 -t 50
#---------------------------------------------------------------------------
echo End Time is `date`
time2=`date +"%Y-%m-%d %H:%M:%S"`
timerun1=$(($(date +%s -d "$time2") - $(date +%s -d "$time1")))
echo $timerun1 | awk '{print "Running time is " $1*1/3600 " hours(" $1*1/60  " mins)"}'
```

Extracting inversion and translocation from .anchor file

```
python3 anchor_trans.sh Oy-0_Yo-0.anchor > Oy-0_Yo-0_inv_trans.txt
```

```
>anchor_trans.sh
import sys,os
anchor = f"{sys.argv[1]}.anchor"
with open(anchor,'r') as f:
	dic_ref = {}
	for x in range(14):
		dic_ref[f"Chr{x+1}"]={}
	for l in f:
		l0=l.strip().split()
		if 'refChr' in l:
			#print(l.strip())
			pass
		else:
			if '#' not in l:
				l0=l.strip().split()
				refChr,referenceStart,referenceEnd,queryChr,queryStart,queryEnd,strand,gene,blockIndex,score = l0[0],l0[1],l0[2],l0[3],l0[4],l0[5],l0[6],l0[7],l0[8],l0[9]
				#if refChr != queryChr:
				if int(referenceEnd) - int(referenceStart) < 0:continue
				pos_key = int(referenceStart)
				dic_ref[refChr][pos_key] = f"{referenceStart}\t{referenceEnd}\t{queryStart}\t{queryEnd}\t{strand}\t{refChr}\t{queryChr}\t{gene}\t{blockIndex}\t{score}"
block = 0
block_dic ={}
for x in dic_ref.keys():
	block_dic[x] = {}
	tmp_region_list_ref,tmp_region_list_alt = [],[]
	for num in range(len(list(dic_ref[x]))):
		if num == len(list(dic_ref[x]))-1:
			list_key_sorted_N0 = sorted(list(dic_ref[x]))[num]
			pos_key_N0 = dic_ref[x][list_key_sorted_N0]
			pos_key_list_N0 = pos_key_N0.split()
			pos_key_strand_N0 = pos_key_list_N0[4]
			ref_chr, alt_chr = pos_key_list_N0[5],pos_key_list_N0[6]
			label_Next = 'Next_Chr_break'  #Break
			tmp_region_list_ref='\t'.join(tmp_region_list_ref)
			tmp_region_list_alt='\t'.join(tmp_region_list_alt)
			block +=1
			if pos_key_strand_N0 == "-":
				inversion = 'Inversion'
			else:
				inversion = 'Normal'
			pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
			print(pri)
			block_dic[x][f"Block{block}"] = pri
			#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,label_Next,pos_key_strand_N0)
		#	print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
			tmp_region_list_ref,tmp_region_list_alt = [],[]
		else:
			list_key_sorted_N0,list_key_sorted_N1 = sorted(list(dic_ref[x]))[num],sorted(list(dic_ref[x]))[num+1]
			pos_key_N0, pos_key_N1 = dic_ref[x][list_key_sorted_N0], dic_ref[x][list_key_sorted_N1]

			pos_key_list_N0, pos_key_list_N1 = pos_key_N0.split(),pos_key_N1.split()
			pos_key_strand_N0, pos_key_strand_N1 = pos_key_list_N0[4], pos_key_list_N1[4]
			ref_chr, alt_chr = pos_key_list_N0[5],pos_key_list_N0[6]
			if pos_key_strand_N0 == pos_key_strand_N1: #strand consistent or not
				if int(pos_key_list_N0[1])+1 == int(pos_key_list_N1[0]): #ref consistent or not
					if pos_key_strand_N0 == '+':
						#if int(pos_key_list_N0[3])+1 == int(pos_key_list_N1[2]): #query consistent or not
						label_Next = 'Next_continue'
						tmp_region_list_ref.append(pos_key_list_N0[0])
						tmp_region_list_ref.append(pos_key_list_N1[1])
						tmp_region_list_ref = [tmp_region_list_ref[0],tmp_region_list_ref[-1]]
						tmp_region_list_alt.append(pos_key_list_N0[2])
						tmp_region_list_alt.append(pos_key_list_N1[3])	
						tmp_region_list_alt = [tmp_region_list_alt[0],tmp_region_list_alt[-1]]
					#		print(num,x,pos_key_N0,'Next_continue',tmp_region_list_ref,tmp_region_list_alt)
						#else:
						#	label_Next = 'Warning1:Next_continueButQueryBreak'
						#	print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
						#	tmp_region_list_ref,tmp_region_list_alt = [],[]
					elif pos_key_strand_N0 == '-':
						#if int(pos_key_list_N1[3])+1 == int(pos_key_list_N0[2]):
						label_Next = 'Next_continue'
						tmp_region_list_ref.append(pos_key_list_N0[0])
						tmp_region_list_ref.append(pos_key_list_N1[1])
						tmp_region_list_ref = [tmp_region_list_ref[0],tmp_region_list_ref[-1]]
						tmp_region_list_alt.append(pos_key_list_N0[3])
						tmp_region_list_alt.append(pos_key_list_N1[2])	
						tmp_region_list_alt = [tmp_region_list_alt[0],tmp_region_list_alt[-1]]
					#		print(num,x,pos_key_N0,'Next_continue',tmp_region_list_ref,tmp_region_list_alt)
						#else:
						#	label_Next = 'Warning2:Next_continueBurQueryBreak'
							#print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
						#	tmp_region_list_ref,tmp_region_list_alt = [],[]
				else:
					label_Next = 'Next_break'  #Break 
					tmp_region_list_ref='\t'.join(tmp_region_list_ref)
					tmp_region_list_alt='\t'.join(tmp_region_list_alt)
					block += 1
					if pos_key_strand_N0 == "-":
						inversion = 'Inversion'
					else:
						inversion = 'Normal'
					pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
					print(pri)
					block_dic[x][f"Block{block}"] = pri
					#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,pos_key_strand_N0,label_Next)
					#print(num,x,pos_key_N0,'Next_break',tmp_region_list_ref,tmp_region_list_alt)
					tmp_region_list_ref,tmp_region_list_alt = [],[]
			elif pos_key_strand_N0 != pos_key_strand_N1:
					label_Next = 'Next_Strand_break'  #Break
					tmp_region_list_ref='\t'.join(tmp_region_list_ref)
					tmp_region_list_alt='\t'.join(tmp_region_list_alt)
					block += 1
					if pos_key_strand_N0 == "-":
						inversion = 'Inversion'
					else:
						inversion = 'Normal'
					pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
					print(pri)
					block_dic[x][f"Block{block}"] = pri
					#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,pos_key_strand_N0,label_Next)
					#print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
					tmp_region_list_ref,tmp_region_list_alt = [],[]
dic_ref={}

print('############################################')
for ch in block_dic.keys():
	if list(block_dic[ch]) == []:continue
	block_list = list(block_dic[ch])
	for i in range(len(block_list)):
		if len(block_list) == 1:
			block0 = block_list[i]
			block_info0 = block_dic[ch][block0]
			block_info_list0 = block_info0.split()
			print(block_info0)
		else:
			if i == 0:
				block0 = block_list[i]
				block1 = block_list[i+1]
				block_info0, block_info1 = block_dic[ch][block0], block_dic[ch][block1]
				block_info_list0, block_info_list1 = block_info0.split(), block_info1.split()
				block0_Qsta, block0_Qend, block1_Qsta, block1_Qend=block_info_list0[5],block_info_list0[6],block_info_list1[5],block_info_list1[6]
				block0_strand, block1_strand = block_info_list0[7], block_info_list1[7]
				ref_chr, alt_chr = block_info_list0[1],block_info_list0[4]
				if ref_chr != alt_chr:
					if block0_strand == "+":
						Trans = "TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
					elif block0_strand == "-":
						Trans = "Inversion/TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
				else:
					if block0_strand == "+":
						if block1_strand == "+":
							judge = int(block1_Qsta) - int(block0_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block1_Qend) - int(block0_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					if block0_strand == "-":
						if block1_strand == "+":
							judge =  int(block1_Qsta) - int(block0_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block1_Qend) - int(block0_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					pass
			else:
				block0 = block_list[i]
				block1 = block_list[i-1]
				block_info0, block_info1 = block_dic[ch][block0], block_dic[ch][block1]
				block_info_list0, block_info_list1 = block_info0.split(), block_info1.split()
				block0_Qsta, block0_Qend, block1_Qsta, block1_Qend=block_info_list0[5],block_info_list0[6],block_info_list1[5],block_info_list1[6]
				block0_strand, block1_strand = block_info_list0[7], block_info_list1[7]
				ref_chr, alt_chr = block_info_list0[1],block_info_list0[4]
				if ref_chr != alt_chr:
					if block0_strand == "+":
						Trans = "TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
					elif block0_strand == "-":
						Trans = "Inversion/TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
				else:

					if block0_strand == "+":
						if block1_strand == "+":
							judge = int(block0_Qsta) - int(block1_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block0_Qsta) - int(block1_Qsta)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					if block0_strand == "-":
						if block1_strand == "+":
							judge =  int(block0_Qend) - int(block1_Qend)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block0_Qend) - int(block1_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)					
```


