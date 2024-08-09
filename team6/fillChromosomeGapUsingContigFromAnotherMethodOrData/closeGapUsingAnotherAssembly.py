import Fasta
from Anchor import readAnchorFile, Anchor
from functools import cmp_to_key
from Gap import compareGapToBeFilled, findAssemblyGaps, GapToBeFilled
import argparse

def closeGapUsingAnotherAssembly(refChrFasta, contigFasta, anchorsFile, outputFile, lineWidth=60):

    refChr_chromosome_names, refChr_fastas = Fasta.readFastaFile(refChrFasta)
    contig_chromosome_names, contig_fastas = Fasta.readFastaFile(contigFasta)

    synBlocks = readAnchorFile(anchorsFile)
    gaps = findAssemblyGaps(refChr_fastas)
    gapsToBeFilled = []
    for gap in gaps:
        for synBlock in synBlocks :
            if gap.chr == synBlock.refChr and gap.start > synBlock.refStart and gap.end < synBlock.refEnd:
                firstAnchor = Anchor
                lastAnchor = Anchor
                notLastAnchorEd = True
                for anchor in synBlock.anchors :
                    if anchor.refEnd < gap.start:
                        firstAnchor = anchor
                    if notLastAnchorEd and anchor.refStart > gap.end:
                        notLastAnchorEd = False
                        lastAnchor = anchor

                gapToBeFilled = GapToBeFilled(firstAnchor, lastAnchor)
                gapsToBeFilled.append(gapToBeFilled)
    gapsToBeFilled = sorted(gapsToBeFilled, key=cmp_to_key(compareGapToBeFilled))
    changed = True
    while changed:
        changed = False
        toberemoved = []
        i = 0
        while i < len(gapsToBeFilled) and (not changed):
            j = i + 1
            while j < len(gapsToBeFilled) and (not changed):
                if (gapsToBeFilled[j].startAnchor.refStart <= gapsToBeFilled[i].startAnchor.refStart and  gapsToBeFilled[i].startAnchor.refStart <= gapsToBeFilled[j].endAnchor.refEnd):
                    changed = True
                    firstAnchor = gapsToBeFilled[j].startAnchor
                    lastAnchor = gapsToBeFilled[j].endAnchor
                    if gapsToBeFilled[i].endAnchor.refEnd > lastAnchor.refEnd:
                        lastAnchor = gapsToBeFilled[i].endAnchor

                    gapToBeFilled = GapToBeFilled(firstAnchor, lastAnchor)
                    gapsToBeFilled.append(gapToBeFilled)
                    toberemoved.append(j)
                    toberemoved.append(i)
                elif (gapsToBeFilled[i].startAnchor.refStart <= gapsToBeFilled[j].startAnchor.refStart and  gapsToBeFilled[j].startAnchor.refStart <= gapsToBeFilled[i].endAnchor.refEnd):
                    changed = True
                    firstAnchor = gapsToBeFilled[i].startAnchor
                    lastAnchor = gapsToBeFilled[j].endAnchor
                    if gapsToBeFilled[i].endAnchor.refEnd > lastAnchor.refEnd:
                        lastAnchor = gapsToBeFilled[i].endAnchor
                    gapToBeFilled = GapToBeFilled(firstAnchor, lastAnchor)
                    gapsToBeFilled.append(gapToBeFilled)
                    toberemoved.append(j)
                    toberemoved.append(i)
                j = j + 1
            i = i + 1
        for index in toberemoved:
            del gapsToBeFilled[index]
    gapsToBeFilled = sorted(gapsToBeFilled, key=cmp_to_key(compareGapToBeFilled), reverse=True)
    for gapToBeFilled in gapsToBeFilled:
        if "+" == gapToBeFilled.startAnchor.strand:
            sequence = Fasta.getSubSequence(refChr_fastas, gapToBeFilled.startAnchor.refChr, 1, gapToBeFilled.startAnchor.refEnd, "+") + Fasta.getSubSequence(contig_fastas, gapToBeFilled.startAnchor.queryChr, gapToBeFilled.startAnchor.queryEnd + 1, gapToBeFilled.endAnchor.queryStart - 1, "+") + Fasta.getSubSequence(refChr_fastas, gapToBeFilled.startAnchor.refChr, gapToBeFilled.endAnchor.refStart, len(refChr_fastas[gapToBeFilled.startAnchor.refChr]), "+")
            refChr_fastas[gapToBeFilled.startAnchor.refChr] = sequence
        elif "-" == gapToBeFilled.startAnchor.strand:
            sequence = Fasta.getSubSequence(refChr_fastas, gapToBeFilled.startAnchor.refChr, 1, gapToBeFilled.startAnchor.refEnd, "+") + Fasta.getSubSequence(contig_fastas, gapToBeFilled.startAnchor.queryChr, gapToBeFilled.endAnchor.queryEnd + 1, gapToBeFilled.startAnchor.queryStart - 1, "-") + Fasta.getSubSequence(refChr_fastas, gapToBeFilled.startAnchor.refChr, gapToBeFilled.endAnchor.refStart, len(refChr_fastas[gapToBeFilled.startAnchor.refChr]), "+")
            refChr_fastas[gapToBeFilled.startAnchor.refChr] = sequence

        print(gapToBeFilled.startAnchor.queryChr + "\t" + str(gapToBeFilled.startAnchor.queryStart) + "\t" +
              str(gapToBeFilled.startAnchor.queryEnd) + "\t" + gapToBeFilled.endAnchor.queryChr + "\t" +
              str(gapToBeFilled.endAnchor.queryStart) + "\t" + str(gapToBeFilled.endAnchor.queryEnd))

        print (gapToBeFilled.startAnchor.refChr + "\t" + str(gapToBeFilled.startAnchor.refStart) + "\t" +
               str(gapToBeFilled.startAnchor.refEnd) + "\t" + gapToBeFilled.endAnchor.refChr + "\t" +
               str(gapToBeFilled.endAnchor.refStart) + "\t" + str(gapToBeFilled.endAnchor.refEnd))

    Fasta.writeFasta(refChr_fastas, outputFile, lineWidth=lineWidth)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Close chromosome level assembly gap using the contig assembly from another assembly software or data')
    parser.add_argument('-r', '--ref', dest='ref_fasta_file', type=str, help="The reference genome of chromosome assembly in FASTA format", required=True)
    parser.add_argument('-q', '--query', dest='query_fasta_file', type=str, help="The query genome of contig or scarfold or chromosome assembly in FASTA format", required=True)
    parser.add_argument('-a', '--anchor', dest='anchor_file', type=str, help="The anchor file generated using AnchorWave", required=True)
    parser.add_argument('-o', '--output', dest='out_file', type=str, help="The output File in FASTA format", required=True)

    args = parser.parse_args()

    closeGapUsingAnotherAssembly(args.ref_fasta_file, args.query_fasta_file, args.anchor_file, args.out_file)



# Example command python3 ../assembly/closeGapUsingAnotherAssembly.py -r chr.7a.fa -q nd.asm.fasta -a align1.anchors -o chr.7a.filled.fa
#
#
# synBlocks = readAnchorFile("/Users/baoxingsong/fillGapUsingContig/07.7a/align1.anchors")
# synBlock = synBlocks[2]
#
# index = 0
# gapsToBeFilled = []
# while index < 10:
#     ii = random.randint(0, len(synBlock.anchors))
#     jj = random.randint(0, len(synBlock.anchors))
#     time.sleep(0.1)
#     if jj > ii:
#         firstAnchor = synBlock.anchors[ii]
#         lastAnchor = synBlock.anchors[jj]
#         gapToBeFilled = GapToBeFilled(firstAnchor, lastAnchor)
#         gapsToBeFilled.append(gapToBeFilled)
#         index = index + 1
#
# print("before sorting")
# for gapToBeFilled in gapsToBeFilled:
#     print(gapToBeFilled.startAnchor.refChr + "\t" + str(gapToBeFilled.startAnchor.refStart) + "\t" + str(gapToBeFilled.startAnchor.refEnd))
#
# gapsToBeFilled = sorted(gapsToBeFilled, key=cmp_to_key(compareGapToBeFilled), reverse=True)
#
# print("")
# print("after sorting")
# for gapToBeFilled in gapsToBeFilled:
#     print(gapToBeFilled.startAnchor.refChr + "\t" + str(gapToBeFilled.startAnchor.refStart) + "\t" + str(gapToBeFilled.startAnchor.refEnd))
