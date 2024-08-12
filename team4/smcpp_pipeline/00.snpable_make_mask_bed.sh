#!/bin/sh
#!/bin/bash
 # usage: ./00.prar_vcffile.sh -s /new/path/to/snpable -p /new/prefix -k 35 -g /new/genome.fa -o /new/output/dir
# 默认值
snpable_script_path="/home/cuixiwang/my_data/biosoft/seqbility-20091110"
prefix="/home/cuixiwang/my_data/hacton/apple.chr"
k=35
GENOME="/home/cuixiwang/my_data/hacton/apple.chr.fa"
OUTDIR="/home/cuixiwang/my_data/hacton"

# 解析命令行参数
while getopts "s:p:k:g:o:" opt; do
  case $opt in
    s) snpable_script_path="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    k) k="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    \?) echo "无效的选项: -$OPTARG" >&2; exit 1 ;;
  esac
done

# 更新PATH
PATH=$PATH:$snpable_script_path

# 打印变量以供调试
echo "snpable_script_path: $snpable_script_path"
echo "prefix: $prefix"
echo "k: $k"
echo "GENOME: $GENOME"
echo "OUTDIR: $OUTDIR"
                    #文件输出位置
mkdir ${OUTDIR}/snpable
cd ${OUTDIR}/snpable
echo "Starting extraction of overlapping ${k}-mer subsequences"
splitfa $GENOME $k | split -l 20000000
cat x* >> ${prefix}_split.$k
# if it can't find splitfa, try adding seqbility to the path using 'PATH=$PATH:/project/WagnerLab/jrick/msmc_Sept2017/snpable/scripts'
echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"
# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
if [ -f "${GENOME}.bwt" ]; then
       echo "$GENOME already indexed"
else
       echo "indexing $GENOME"
       bwa index $GENOME
fi
echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 66 -R 1000000 -O 3 -E 3 ${GENOME} ${prefix}_split.${k} > ${prefix}_split.${k}.sai
bwa samse -f ${prefix}_split.${k}.sam $GENOME ${prefix}_split.${k}.sai ${prefix}_split.${k}
echo "reads aligned, starting to generate rawMask"
gen_raw_mask.pl ${prefix}_split.${k}.sam > ${prefix}_rawMask.${k}.fa
echo "raw mask created as ${prefix}_rawMask.35.fa, now generating final mask with stringency r=60%"
gen_mask -l ${k} -r 0.6 ${prefix}_rawMask.${k}.fa > ${prefix}_mask.${k}.60.fa
echo "all done! final mask saved as ${prefix}_mask.${k}.60.fa"