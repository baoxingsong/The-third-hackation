#!/bin/python3

# script from: https://raw.githubusercontent.com/twooldridge/misc/master/SMC_bootstrap_BW.py

# This script is a modification of the SUPER HELPFUL script found in this github issue: https://github.com/popgenmethods/smcpp/issues/37.
# I've modified it to deal with irregularities in the input file, say due to some of the features of masking that vcf2smc produces
import click
import sys
import subprocess
import os
import random
import gzip
import time

@click.command()
@click.option('--nr_bs', type=int, help="nr of bs [20]", default=20)
@click.option("--chunk_size", type=int, help="size of bootstrap chunks [5000000]", default=5000000)
@click.option("--chunks_per_chromosome", type=int, help="nr of chunks to put on one chromosome in the bootstrap [20]", default=20)
@click.option("--nr_chromosomes", type=int, help="nr of chromosomes to write [30]", default=24)
@click.option("--seed", type=int, help="initialize the random number generator", default=None)
@click.argument("out_dir_prefix")
@click.argument("files", nargs=-1)
def bs(nr_bs, chunk_size, chunks_per_chromosome, nr_chromosomes, seed, out_dir_prefix, files):
    chunks = []
    offset = 0
    chunks_in_chrom = []
    if not seed:
        seed = int(time.time())
    random.seed(seed)
    print('seed: %s' % seed)
    for file in files:
        print(file)
        with gzip.open(file, 'rb') as f:
            header = f.readline().decode()
            prev_chunk_idx = -1
            for line_count, line in enumerate(f):
                line = line.decode()
                line = [int(x) for x in line.strip().split()]
                if line_count == 0:
                    pos = 1
                    offset += 1
                else:
                    pos = line[0] + offset
                    offset += line[0]
                chunk_index = (pos - 1) // chunk_size
                if (chunk_index - prev_chunk_idx) >= 2:
                    print('\n')
                    print('Next position is in %s bases' % pos)
                    print('Last relative position was at %s bases' % lastpos)
                    print('This might be the result of how vcf2smc masks large stretches of the vcf (say for subsampling a chromosome), or it might be a real error.\nDouble check that this input file is formatted how you intend.\nStopping "chunking" of input file and moving on\n')
                    break
                print('Processing chunk %s of input file' % chunk_index)
                if chunk_index > len(chunks_in_chrom) - 1:
                    chunks_in_chrom.append([])
                chunks_in_chrom[chunk_index].append(line)
                prev_chunk_idx = prev_chunk_idx + (chunk_index - prev_chunk_idx)
                lastpos = pos

        chunks.extend(chunks_in_chrom)

    for bootstrap_id in range(1, nr_bs + 1):
        for chr_ in range(1, nr_chromosomes + 1):
            chr_dir = "{}_{}".format(out_dir_prefix, bootstrap_id)
            if not os.path.exists(chr_dir):
                os.makedirs(chr_dir)
            chr_file = "{}/bootstrap_chr{}.gz".format(chr_dir, chr_)
            print("writing", chr_file, file=sys.stderr)
            with gzip.open(chr_file, 'wb') as f:
                f.write(header.encode())
                for i in range(chunks_per_chromosome):
                    chunk_id = random.randrange(len(chunks))
                    for line in chunks[chunk_id]:
                        line = ' '.join([str(x) for x in line]) + '\n'
                        f.write(line.encode())

@click.command()
def dg():
    if len(sys.argv) != 5:
        print("Usage: python script.py dg vcf_file POPName ID_file distinguish_file")
        sys.exit(1)

    vcf_file = sys.argv[2]
    POP = sys.argv[3]
    ID_file = sys.argv[4]
    dis = sys.argv[5]

    # 读取样本ID并用逗号分隔
    with open(ID_file, 'r') as f:
        line = ','.join([x.strip() for x in f.readlines()])
    print(line)

    # 打开dged_file
    with open(dis, 'r') as f:
        dged_lines = [x.strip() for x in f.readlines()]

    for d in dged_lines:
        print(d)
        i = 1
        while i < 10:
            # 执行smc++ vcf2smc命令
            cmd = f"smc++ vcf2smc  -d {d} {d} {vcf_file} {d}.Chr{i}.smc.gz Chr{i} {POP}:{line}"
            subprocess.run(cmd, shell=True)
            # 输出生成的SMC文件名
            print(f"{d}.Chr{i}.smc.gz")
            i += 1

@click.group()
def cli():
    pass

cli.add_command(bs)
cli.add_command(dg)

if __name__ == "__main__":
    cli()
