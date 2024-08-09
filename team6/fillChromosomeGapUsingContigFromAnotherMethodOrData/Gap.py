import Anchor
class Gap:
    chr = ""
    start = 0
    end = 0
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

def findAssemblyGaps (fastas):
    gaps = []
    for chr in fastas:
        prev = ' '
        i = 0
        start = 0
        while i < len(fastas[chr]):
            if fastas[chr][i] == 'N':
                if prev != 'N':
                    start = i + 1 # the is 0 based position, and start is 1 based position
            else:
                if prev == 'N':
                    gap = Gap(chr, start, i) # end = i ,  end is 1 based position, and both and end are included
                    gaps.append(gap)

            prev = fastas[chr][i]
            i = i + 1

        if prev == 'N': # if the end with 'N'
            gap = Gap(chr, start, i)  # end = i ,  end is 1 based position, and both and end are included
            gaps.append(gap)

    return gaps

class GapToBeFilled:
    startAnchor = Anchor
    endAnchor = Anchor
    def __init__(self, startAnchor, endAnchor):
        self.startAnchor = startAnchor
        self.endAnchor = endAnchor

def compareGapToBeFilled(gapToBeFilled1, gapToBeFilled2):
    if gapToBeFilled1.startAnchor.refChr == gapToBeFilled2.startAnchor.refChr:
        if gapToBeFilled1.startAnchor.refStart < gapToBeFilled2.startAnchor.refStart:
            return -1
        else:
            return 1
    elif gapToBeFilled1.startAnchor.refChr < gapToBeFilled2.startAnchor.refChr:
        return -1
    elif gapToBeFilled1.startAnchor.refChr > gapToBeFilled2.startAnchor.refChr:
        return 1
    else:
        return 0



def compareGaps(gap1, gap2):
    if gap1.chr == gap2.chr:
        if gap1.start < gap2.start:
            return -1
        else:
            return 1
    elif gap1.chr < gap2.chr:
        return -1
    elif gap1.chr > gap2.chr:
        return 1
    else:
        return 0