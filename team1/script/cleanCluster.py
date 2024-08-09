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

def giveGeneAnId(gffFile):
    chromosome_gene_dict, chromosome_gene_list, geneName_toChr_dict, transcript_gene_map = readGff(gffFile)
    geneIndex = dict()

    for chr in chromosome_gene_list:
        i = 0
        while i < len(chromosome_gene_list[chr]):
            geneIndex[chromosome_gene_list[chr][i].name] = i + 1
            i += 1

    return geneIndex, geneName_toChr_dict, transcript_gene_map

def readCluster(clusterFile):
    cluster = dict()
    geneNameWithPreFixToGeneName = dict()
    with open(clusterFile) as f:
        for line in f:
            elements = line.split()
            geneId = elements[0]
            clusterId = elements[1]
            if clusterId not in cluster:
                cluster[clusterId] = dict()
            m = re.search(r'^(\S+?)_(\S+)', line)
            if m != None:
                species = m.group(1)
                geneNameWithPreFixToGeneName[geneId] = m.group(2)
                if species not in cluster[clusterId]:
                    cluster[clusterId][species] = dict()
                    cluster[clusterId][species] = []

                cluster[clusterId][species].append(geneId)
    return cluster, geneNameWithPreFixToGeneName


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prepare input file for AnchorWave collinearity gene analysis')
    parser.add_argument('-c', '--clusterFile', dest='clusterFile', type=str, help="clusterFile file", required=True)
    parser.add_argument('-m', '--maximumDistance', dest='maximumDistance', type=int, default=5,  help="maximumDistance allowed as tandem, longer than this distance, the whole cluster will be dropped")
    parser.add_argument('-g', '--gffFiles', dest='gffFiles', type=str, help="prefix1:gffFilePath1;prefix2:gffFilePath2;", required=True)
    parser.add_argument('-o', '--output', dest='out_file', type=str, help="output File", required=True)

    args = parser.parse_args()
    cluster, geneNameWithPreFixToGeneName = readCluster(args.clusterFile)

    genesIndex = dict()
    geneName_toChr_dicts = dict()
    transcript_gene_maps = dict()
    elements = args.gffFiles.split(';')
    for element in elements:
        es = element.split(':')
        geneIndex, geneName_toChr_dict, transcript_gene_map = giveGeneAnId(es[1])
        genesIndex[es[0]] = geneIndex
        geneName_toChr_dicts[es[0]] = geneName_toChr_dict
        transcript_gene_maps[es[0]] = transcript_gene_map

    target_output = open(args.out_file, 'w')
    for clusterid in cluster:
        useThisCluster = True
        for species in cluster[clusterid]:
            if len(cluster[clusterid][species]) > 1:
                ids = []
                chrs = set()
                for gene in cluster[clusterid][species]:
                    ids.append(genesIndex[species][transcript_gene_maps[species][geneNameWithPreFixToGeneName[gene]]])
                    chrs.add(geneName_toChr_dicts[species][transcript_gene_maps[species][geneNameWithPreFixToGeneName[gene]]])
                if len(chrs) > 1:
                    useThisCluster = False
                ids = sorted(set(ids))
                idIndex = 0
                while idIndex < len(ids) - 1:
                    idIndex2 = idIndex + 1
                    if ids[idIndex2] - ids[idIndex] > args.maximumDistance:
                        useThisCluster = False
                    idIndex = idIndex + 1

        if useThisCluster :
            for species in cluster[clusterid]:
                if len(cluster[clusterid][species]) > 0:
                    target_output.write( cluster[clusterid][species][0] + "\t" + clusterid + "\n")
    target_output.close()
