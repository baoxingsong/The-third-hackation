2. ABBABABAwindows.py
git clone https://github.com/simonhmartin/genomics_general   #download the conmands from the github
python parseVCF.py -i SNP.vcf.gz --skipIndels --minQual 30  |bgzip > output.geno.gz # deal the format of the SNP file


