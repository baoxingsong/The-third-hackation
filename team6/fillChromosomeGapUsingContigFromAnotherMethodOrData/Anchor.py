import this


class Anchor:
    refChr = ""
    refStart = 0
    refEnd = 0
    queryChr = ""
    queryStart = 0
    queryEnd = 0
    strand = ""

    def __init__(self, refChr, refStart, refEnd, queryChr, queryStart, queryEnd, strand):
        self.refChr = refChr
        self.refStart = refStart
        self.refEnd = refEnd
        self.queryChr = queryChr
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.strand = strand
class Block:
    anchors = []
    refStart = 0
    refEnd = 0
    queryStart = 0
    queryEnd = 0
    refChr = ""
    queryChr = ""
    def __init__(self):
        self.anchors = []

    def updateInfor(self):
        self.refChr = self.anchors[0].refChr
        self.queryChr = self.anchors[0].queryChr

        self.refStart = self.anchors[0].refStart
        if self.refStart > self.anchors[0].refEnd:
            self.refStart = self.anchors[0].refEnd

        if self.refStart > self.anchors[ len(self.anchors)-1 ].refEnd:
            self.refStart = self.anchors[ len(self.anchors)-1].refEnd

        if self.refStart > self.anchors[len(self.anchors)-1].refStart:
            self.refStart = self.anchors[len(self.anchors)-1].refStart


        self.refEnd = self.anchors[0].refStart
        if self.refEnd < self.anchors[0].refEnd:
            self.refEnd = self.anchors[0].refEnd

        if self.refEnd < self.anchors[len(self.anchors) - 1].refEnd:
            self.refEnd = self.anchors[len(self.anchors) - 1].refEnd

        if self.refEnd < self.anchors[len(self.anchors) - 1].refStart:
            self.refEnd = self.anchors[len(self.anchors) - 1].refStart


        self.queryStart = self.anchors[0].queryStart
        if self.queryStart > self.anchors[0].queryEnd:
            self.queryStart = self.anchors[0].queryEnd

        if self.queryStart > self.anchors[len(self.anchors) - 1].queryEnd:
            self.queryStart = self.anchors[len(self.anchors) - 1].queryEnd

        if self.queryStart > self.anchors[len(self.anchors) - 1].queryStart:
            self.queryStart = self.anchors[len(self.anchors) - 1].queryStart

        self.queryEnd = self.anchors[0].queryStart
        if self.queryEnd < self.anchors[0].queryEnd:
            self.queryEnd = self.anchors[0].queryEnd

        if self.queryEnd < self.anchors[len(self.anchors) - 1].queryEnd:
            self.queryEnd = self.anchors[len(self.anchors) - 1].queryEnd

        if self.queryEnd < self.anchors[len(self.anchors) - 1].queryStart:
            self.queryEnd = self.anchors[len(self.anchors) - 1].queryStart

def readAnchorFile(anchorFile):
    blocks = []
    block = Block()
    with open(anchorFile) as f:
        for line in f:
            line = line.strip()
            elements = line.split()
            if line.startswith("#block begin"):
                block = Block()

            if line.startswith("#block end"):
                block.updateInfor()
                blocks.append(block)

            if line[0] != "#" and line[0:6] != "refChr" and elements[7] != "interanchor" and elements[7][0:14] != "localAlignment":
                anchor = Anchor(elements[0], int(elements[1]), int(elements[2]), elements[3], int(elements[4]), int(elements[5]), elements[6])
                block.anchors.append(anchor)
    return blocks

#synBlocks = readAnchorFile("/Users/baoxingsong/fillGapUsingContig/07.7a/align1.anchors")

#for snynBlock in synBlocks:
#    print(snynBlock.refChr + "\t" + str(snynBlock.refStart) + "\t" + str(snynBlock.refEnd) + "\t" + snynBlock.queryChr + "\t" + str(snynBlock.queryStart) + "\t" + str(snynBlock.queryEnd))
