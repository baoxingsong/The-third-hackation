import re

import Fasta
from functools import cmp_to_key
from Gap import Gap, findAssemblyGaps, compareGaps
import argparse


def deleteGapUpAndDownStream (inputFile, outFile, extendBps = 500000, replace=False):
    chromosome_names, fastas = Fasta.readFastaFile(inputFile)
    gaps = findAssemblyGaps(fastas)
    newGaps = []
    largestGap  = 0
    for gap in gaps:
        print("oldGap\t" + gap.chr + "\t" + str(gap.start) + "\t" + str(gap.end))
        thisGapSize = gap.end - gap.start + 1
        if thisGapSize > largestGap:
            largestGap = thisGapSize
        newStart = gap.start-extendBps
        if newStart < 1:
            newStart = 1
        newEnd = gap.end+extendBps
        if newEnd > len(fastas[gap.chr]):
            newEnd = len(fastas[gap.chr])
        newGap = Gap(gap.chr, newStart, newEnd)
        newGaps.append(newGap)

    newGaps = sorted(newGaps, key=cmp_to_key(compareGaps), reverse=True)

    for gap in newGaps:
        print("newGap\t" + gap.chr + "\t" + str(gap.start) + "\t" + str(gap.end))
        length = gap.end - gap.start + 1
        sequence = Fasta.getSubSequence(fastas, gap.chr, 1, gap.start - 1, "+") + "N" * length + Fasta.getSubSequence(fastas, gap.chr, gap.end + 1, len(fastas[gap.chr]), "+")
        fastas[gap.chr] = sequence

    if replace:
        for chr in fastas:
            print( chr + "\told length:" + str(len(fastas[chr])))
            sequence = fastas[chr]
            pattern = "N{" + str(largestGap) + ",}"
            sequence = re.sub(pattern, "N"*largestGap, sequence)
            fastas[chr] = sequence
            print("new length:" + str(len(fastas[chr])))
            print()

    Fasta.writeFasta(fastas, outFile, lineWidth=60)

# python3 ../assembly/deleteGapUpAndDownStream.py -i chr.7a.fa -e 500000 -o chr.7a.trimed.gap.fa
# python3 ../assembly/deleteGapUpAndDownStream.py -i chr.7a.fa -e 500000 -o chr.7a.trimed.removed.gap.fa -r True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mask the neaby regions of genome assembly gaps.')
    parser.add_argument('-i', '--input', dest='input', type=str, help="The input file with gaps marked as Ns in FASTA format", required=True)
    parser.add_argument('-o', '--output', dest='output', type=str, help="The output File in FASTA format", required=True)
    parser.add_argument('-e', '--extend', dest='extend', type=int, default=500000, help="The number of base pairs to be extended. The extended basepairs will be maksed as 'Ns'")
    parser.add_argument('-r', '--replace', dest='replace', type=bool, default=False, help="Whether trim the masked region (trim as the X 'N's, X is the size of the largest gap) ")
    args = parser.parse_args()
    deleteGapUpAndDownStream(args.input,  args.output, args.extend, args.replace)
