#!/bin/sh
# 默认值
chr=18
use_mask=false
boots=20

# 处理命令行选项
while getopts "chr:m:boots" opt; do
  case $opt in
    chr)
      chr=$OPTARG
      ;;
    m)
      use_mask=true
      ;;
	boots)
      boots=$OPTARG
      ;;
    \?)
      echo "无效的选项: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "选项 -$OPTARG 需要一个参数" >&2
      exit 1
      ;;
  esac
done

# 移除已处理的选项
shift $((OPTIND - 1))

if [ "$1" = "--help" ]; then
 echo "Usage: $0 [-chr <number>] [-m] <poplabel> <distinguished> <idlist|sep by \n>"
 echo "Description: This script requires exactly 3 arguments."
 echo "Arguments:"
 echo "  <poplabel>       Description of poplabel"
 echo "  <distinguished>  Description of distinguished"
 echo "  <idlist>        Description of idname, separated by \n"
 echo "Options:"
 echo "  -chr <number>    Define chr as 1 to <number> (default is 29)"
 echo "  -m               Use mask file in smc++ vcf2smc command"
<<<<<<< HEAD
 echo "  -boost <number>  Define bootstraps number (default is 20)"
=======
 echo "  -boost <number>          Define bootstraps number (default is 20)"
>>>>>>> 80260d2910b942a1ca81090b1c4540c1eb14fb9d

 exit 0
fi


export poplabel=$1
export distinguished=$2
export idlist=$3 

idname=$(echo "$idlist" | tr '\n' ',' | sed 's/,$//')  # list>>>seq ,<idname|sep by commas(,)>
#converts VCF data to the format used by SMC++
# 使用 for 循环从 1 到 chr
for i in $(seq 1 $chr)
do
  if [ "$use_mask" = true ]; then
    smc++ vcf2smc -d $distinguished $distinguished --mask $i.mask.bed.gz $i.vcf.gz $distinguished.$i.smc.gz $i $poplabel:$idname
  else
    smc++ vcf2smc -d $distinguished $distinguished $i.vcf.gz $distinguished.$i.smc.gz $i $poplabel:$idname
  fi
done
fi

if [ "$boots" -gt 0 ]; then
  echo "Use boots" #bootstrapping by breaking up the genome into 5-Mb segments and then randomly sampling with replacement
  bootstrap_smcpp.py bs --nr_chromosomes $chr  --nr_bootstraps $boots  $distinguished $distinguished.*.smc.gz 
  smc++ estimate --cores 30 -o $out -o bs_n/ 3.9e-9 PATH/*.smc.gz --spline cubic
else
  smc++ estimate --cores 30 -o $out -o bs_n/ 3.9e-9 PATH/*.smc.gz --spline cubic
fi