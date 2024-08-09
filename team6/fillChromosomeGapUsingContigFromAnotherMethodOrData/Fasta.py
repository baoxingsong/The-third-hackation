import re

# baoxing.song@pku-iaas.edu.cn

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
                    fastas[name] = s
                    chromosome_names.append(name)
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            fastas[name] = s
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
    if start > len(fastas[name]):
        return ""
    if end > len(fastas[name]):
        end = len(fastas[name])

    seq = fastas[name][start:end]
    if "+" == strand:
        return seq
    else:
        return getReverseComplementary(seq)


def writeFasta(fastas, file, lineWidth = 80):
    f = open(file, "w")
    for fasta in fastas:
        f.write(">" + fasta + "\n")
        start = 0
        while start < len(fastas[fasta]):
            end = start + lineWidth
            if end > len(fastas[fasta]):
                end = len(fastas[fasta])
            f.write(fastas[fasta][start:end] + "\n")
            start = end
    f.close()
