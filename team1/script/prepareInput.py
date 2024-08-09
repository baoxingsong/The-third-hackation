import argparse
import os
import re
import numpy as np
class Transcript:
    chromosome_name = ""
    if_orf_conserved = True
    meta_informaiton = ""
    name = ""
    start = 1000000000
    end = 0
    strand = ""
    genome_sequence = ""
    cds_sequence = ""
    exons = np.empty([0, 2], int)  # start and end position
    Cds = np.empty([0, 2], int)

    def __init__(self, name, strand, chromosome_name):
        self.name = name
        self.strand = strand
        self.chromosome_name = chromosome_name

    def add_cds(self, start, end):
        self.Cds = np.append(self.Cds, np.array([[start, end]]), axis=0)
        # self.Cds.append([start, end])

    def add_exon(self, start, end):
        self.exons = np.append(self.exons, np.array([[start, end]]), axis=0)

    def updateCordinate(self):
        if "+" == self.strand:
            self.Cds = np.sort(self.Cds, axis=0)
        else:
            self.Cds = -np.sort(-self.Cds, axis=0)
        if len(self.Cds) > 0:
            if "+" == self.strand:
                self.start = self.Cds[0][0]
                self.end = self.Cds[len(self.Cds)-1][1]
            else:
                self.start = self.Cds[len(self.Cds) - 1][0]
                self.end = self.Cds[0][1]

    def __lt__(self, other):
        if self.start < other.start:
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if self.start < other.start:
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False


class Gene:
    name = ""
    start = 100000000
    end = 0
    transcripts = np.empty([0, 1], Transcript)  # start and end position
    length = 0
    strand = ""

    def __init__(self, name, strand):
        self.name = name
        self.strand = strand

    def add_transcript(self, transcript):
        self.transcripts = np.append(self.transcripts, [transcript])

    def updateCordinate(self):
        self.transcripts = np.sort(self.transcripts)
        for transcript in self.transcripts:
            self.start = self.transcripts[0].start
            if transcript.end > self.end:
                self.end = transcript.end
        self.length = self.end - self.start

    def __lt__(self, other):
        if self.start < other.start:
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if self.start < other.start:
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False
# end of gene class


def readGff(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_dict = dict()
    geneName_toChr_dict = dict()
    chromosome_transcript_dict = dict()
    chromosome_gene_list = dict()  # start and end position
    fake_transcript_gene_map = dict()
    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if m == None:
                m = re.search(r'^(\S+)\t(\S+)\tCDS\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([:.\-\w]+?);', line)
                if m == None:
                    m = re.search(r'^(\S+)\t(\S+)\tCDS\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([:.\-\w]+?)$', line)
                if m != None:
#                    print(line)
                    chromosome_name = m.group(1)
                    if chromosome_name not in chromosome_transcript_dict:
                        chromosome_transcript_dict[chromosome_name]=dict()

                    start = int(m.group(3))
                    end = int(m.group(4))
                    if start > end:
                        start, end = end, start

                    transcript_name = m.group(8)
                    transcript_name = transcript_name.replace("transcript:", "")
                    if transcript_name not in chromosome_transcript_dict[chromosome_name]:
                        strand = m.group(6)
                        transcript = Transcript(transcript_name, strand, chromosome_name)
#                        print(transcript_name + "\t" + strand + "\t" + chromosome_name)
                        (chromosome_transcript_dict[chromosome_name])[transcript_name] = transcript
                    chromosome_transcript_dict[chromosome_name][transcript_name].add_cds(start, end)

                m = re.search(r'^(\S+)\t(\S+)\t\S+\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([:.\-\w]+?);Parent=([:.\-\w]+?);', line)
                if m == None:
                    m = re.search(r'^(\S+)\t(\S+)\t\S+\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([:.\-\w]+?);Parent=([:.\-\w]+?)$', line)
                if m != None:
                    fake_transcript_name = m.group(8)
                    fake_gene_name = m.group(9)
                    fake_transcript_name = fake_transcript_name.replace("transcript:", "")
                    fake_transcript_gene_map[fake_transcript_name] = fake_gene_name

    # update transcript, get the start end information
    for chromosome_name in chromosome_transcript_dict:
        if chromosome_name not in chromosome_gene_dict:
            chromosome_gene_dict[chromosome_name] = dict()
        for transcript_name in chromosome_transcript_dict[chromosome_name]:
            transcript = chromosome_transcript_dict[chromosome_name][transcript_name]
            transcript.updateCordinate()
            gene_name = fake_transcript_gene_map[transcript_name]
            if gene_name not in chromosome_gene_dict[chromosome_name]:
                chromosome_gene_dict[chromosome_name][gene_name] = Gene(gene_name, transcript.strand)

            chromosome_gene_dict[chromosome_name][gene_name].add_transcript(transcript)
    # print("transcript update done")
    for chromosome_name in chromosome_gene_dict:
        gene_list = np.empty([0, 1], Gene)
        for gene_name in chromosome_gene_dict[chromosome_name]:
            chromosome_gene_dict[chromosome_name][gene_name].updateCordinate()
            gene_list = np.append(gene_list, [chromosome_gene_dict[chromosome_name][gene_name]])
            geneName_toChr_dict[gene_name] = chromosome_name
        gene_list = np.sort(gene_list)
        chromosome_gene_list[chromosome_name] = gene_list

    return chromosome_gene_dict, chromosome_gene_list, geneName_toChr_dict, fake_transcript_gene_map


def update_sequence_information_onechromosome(fastas, chromosome_gene_dict, chromosome_name):
    if chromosome_name in fastas:
        if chromosome_name in chromosome_gene_dict:
            for gene_name in chromosome_gene_dict[chromosome_name]:
                for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                    transcript.genome_sequence = getSubSequence(fastas, chromosome_name, transcript.start, transcript.end, transcript.strand)
                    cds_seq_list = []
                    for cds in transcript.Cds:
                        cds_seq_list.append(
                            getSubSequence(fastas, chromosome_name, cds[0], cds[1], transcript.strand))
                    transcript.cds_sequence = "".join(cds_seq_list)


def update_sequence_information(fastas, chromosome_gene_dict):
    for chromosome_name in fastas:
        update_sequence_information_onechromosome(fastas, chromosome_gene_dict, chromosome_name)




class Fasta:
    name = ""
    seq = ""

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq


def readFastaFile(fastaFile):
    fastas = {}
    chromosome_names = []
    name = ""
    seq = []
    with open(fastaFile) as f:
        for line in f:
            m = re.search(r'^>(\S+)', line)
            if m != None:
                if (len(name) > 0) & (len(seq) > 0):
                    s = ''.join(seq)
                    s = re.sub("\\s", "", s)
                    s = s.upper()
                    fasta = Fasta(name, s)
                    fastas[name] = fasta
                    chromosome_names.append(name)
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            fasta = Fasta(name, s)
            fastas[name] = fasta
            chromosome_names.append(name)

    return chromosome_names, fastas


def getReverseComplementary(sequence):
    reversecomplementary = []
    for c in sequence[::-1]:
        if 'A' == c:
            c = 'T'
        elif 'T' == c:
            c = 'A'
        elif 'U' == c:
            c = 'A'
        elif 'C' == c:
            c = 'G'
        elif 'G' == c:
            c = 'C'
        elif 'R' == c:
            c = 'Y'
        elif 'Y' == c:
            c = 'R'
        elif 'K' == c:
            c = 'M'
        elif 'M' == c:
            c = 'K'
        elif 'B' == c:
            c = 'V'
        elif 'V' == c:
            c = 'B'
        elif 'D' == c:
            c = 'H'
        elif 'H' == c:
            c = 'D'
        reversecomplementary.append(c)

    return ''.join(reversecomplementary)


def getSubSequence(fastas, name, start, end, strand):
    # get a sequence fragment from fasta records
    start = start - 1
    if start > len(fastas[name].seq):
        return ""
    if end > len(fastas[name].seq):
        end = len(fastas[name].seq)

    seq = fastas[name].seq[start:end]
    if "+" == strand:
        return seq
    else:
        return getReverseComplementary(seq)

def anchorwave_quota(refGffFile, queryGffFile, blastpresult, outputFile, bit_score, align_length,  reference_prefix, query_prefix):
    target_output = open(outputFile, 'w')
    refChromosome_gene_dict, refChromosome_gene_list, ref_GeneName_toChr_dict, ref_transcript_gene_map = readGff(refGffFile)
    queryChromosome_gene_dict, queryChromosome_gene_list, query_GeneName_toChr_dict, qeury_transcript_gene_map = readGff(queryGffFile)

    refGeneIndex = dict()

    for refChr in refChromosome_gene_list:
        i = 0
        while i < len(refChromosome_gene_list[refChr]):
            refGeneIndex[refChromosome_gene_list[refChr][i].name] = i + 1
            i += 1

    queryGeneIndex = dict()

    for queryChr in queryChromosome_gene_list:
        i = 0
        while i < len(queryChromosome_gene_list[queryChr]):
            queryGeneIndex[queryChromosome_gene_list[queryChr][i].name] = i + 1
            i += 1

    match_pairs = set()
    with open(blastpresult) as f:
        for line in f:
            elements = line.split()
            qseqid = elements[0]
            sseqid = elements[1]

            qseqid = qseqid.replace(query_prefix, "")
            sseqid = sseqid.replace(reference_prefix, "")

            pident = str(float(elements[2]))
            length = int(elements[3])
            mismatch = elements[4]
            gapopen = elements[5]
            qstart = elements[6]
            qend = elements[7]
            sstart = elements[8]
            send = elements[9]
            evalue = elements[10]
            bitscore = float(elements[11])

            match_pair = sseqid + "_" + qseqid
            if (match_pair not in match_pairs) and (bitscore > float(bit_score)) and (length > float(align_length)):
                match_pairs.add(match_pair)
                target_output.write(
                    reference_prefix + sseqid + "\t" + ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]] + "\t" + str(refGeneIndex[ref_transcript_gene_map[sseqid]]) + "\t"
                    + str(refChromosome_gene_dict[ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]]][ref_transcript_gene_map[sseqid]].start) + "\t"
                    + str(refChromosome_gene_dict[ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]]][ref_transcript_gene_map[sseqid]].end) + "\t"
                    + refChromosome_gene_dict[ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]]][ref_transcript_gene_map[sseqid]].strand + "\t"
                    + query_prefix  + qseqid + "\t" + query_GeneName_toChr_dict[qeury_transcript_gene_map[qseqid]] + "\t" + str(queryGeneIndex[qeury_transcript_gene_map[qseqid]]) + "\t" +
                    str(queryChromosome_gene_dict[query_GeneName_toChr_dict[qeury_transcript_gene_map[qseqid]]][qeury_transcript_gene_map[qseqid]].start) + "\t" +
                    str(queryChromosome_gene_dict[query_GeneName_toChr_dict[qeury_transcript_gene_map[qseqid]]][qeury_transcript_gene_map[qseqid]].end) + "\t" +
                    queryChromosome_gene_dict[query_GeneName_toChr_dict[qeury_transcript_gene_map[qseqid]]][qeury_transcript_gene_map[qseqid]].strand + "\t" + pident + "\n")

    target_output.close()
#
# refChromosome_gene_dict, refChromosome_gene_list, ref_GeneName_toChr_dict, ref_transcript_gene_map = readGff("/Users/baoxingsong/Zema.gff")
# refGeneIndex = dict()
#
# for refChr in refChromosome_gene_list:
#     i = 0
#     while i < len(refChromosome_gene_list[refChr]):
#         refGeneIndex[refChromosome_gene_list[refChr][i].name] = i + 1
#         i += 1
#
# sseqid = 'Zm00014ba000020_T001'
# #sseqid = 'Zm00014ba000010'
# print( ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]])

#
# queryChromosome_gene_dict, queryChromosome_gene_list, query_GeneName_toChr_dict, qeury_transcript_gene_map = readGff("/Users/baoxingsong/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3")
#
# queryGeneIndex = dict()
#
# for queryChr in queryChromosome_gene_list:
#     i = 0
#     while i < len(queryChromosome_gene_list[queryChr]):
#         queryGeneIndex[queryChromosome_gene_list[queryChr][i].name] = i + 1
#         i += 1
#
# sseqid = 'KQK86056'
# print( ref_GeneName_toChr_dict[ref_transcript_gene_map[sseqid]])
#
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prepare input file for AnchorWave collinearity gene analysis')
    parser.add_argument('-r', '--ref', dest='ref_gff_file', type=str, help="reference genome GFF file", required=True)
    parser.add_argument('-q', '--query', dest='query_gff_file', type=str, help="query genome GFF file", required=True)

    parser.add_argument('-rp',  dest='reference_prefix', type=str, help="prefix of reference", required=True)
    parser.add_argument('-qp',  dest='query_prefix', type=str, help="prefix of query", required=True)

    parser.add_argument('-b', '--blast_result', dest='blast_result',type=str, help="blast result", required=True)
    parser.add_argument('-o', '--output', dest='out_file', type=str, help="output File", required=True)
    parser.add_argument('-s', '--bitscore', dest='bitscore', type=int, default=100, help="minimum bit score that counts")
    parser.add_argument('-l', '--align_length', dest='align_length', type=int, default=100, help="minimum alignment length that count")

    args = parser.parse_args()
    anchorwave_quota(args.ref_gff_file, args.query_gff_file, args.blast_result,
                                                      args.out_file, args.bitscore, args.align_length, args.reference_prefix, args.query_prefix)
