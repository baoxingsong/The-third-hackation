For chromosome level genome assembly tasks, the contigs/scarfolds are generally linked with Ns.

# deleteGapUpAndDownStream.py
The sequences nearby the gaps are generally not accurant, the deleteGapUpAndDownStream.py is designed to mask or delete the sequences nearby the gaps.
example commands are:

```
python3 ../assembly/deleteGapUpAndDownStream.py -i chr.7a.fa -e 500000 -o chr.7a.trimed.gap.fa
python3 ../assembly/deleteGapUpAndDownStream.py -i chr.7a.fa -e 500000 -o chr.7a.trimed.removed.gap.fa -r True
```
The description about the parameters are :

```
usage: deleteGapUpAndDownStream.py [-h] -i INPUT -o OUTPUT [-e EXTEND] [-r REPLACE]

Mask the neaby regions of genome assembly gaps.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input file with gaps marked as Ns in FASTA format
  -o OUTPUT, --output OUTPUT
                        The output File in FASTA format
  -e EXTEND, --extend EXTEND
                        The number of base pairs to be extended. The extended basepairs will be maksed as 'Ns'
  -r REPLACE, --replace REPLACE
                        Whether trim the masked region (trim as the X 'N's, X is the size of the largest gap)
```

# closeGapUsingAnotherAssembly.py
Then we could use the AnchorWave software (the proali command) to aligned the contigs from another methods or method to the chromosome assembly, and fill the gap using the sequences from the contig assembly.

For example:
`chr.7a.fa` is the chromosome level assembly using Hifi reads and HIC data. \
`zz.t2.gff3` is the annotation for chromosome level assembly. \
`nd.asm.fasta` is the contig assembly using ONT reads.

```
anchorwave gff2seq -r chr.7a.fa -i zz.t2.gff3 -o cds.fa
minimap2 -x splice -t 128 -k 12 -a -p 0.4 -N 20 chr.7a.fa cds.fa > ref.sam
minimap2 -x splice -t 128 -k 12 -a -p 0.4 -N 20 nd.asm.fasta cds.fa > cds.sam
anchorwave proali -i zz.t2.gff3 -r chr.7a.fa -a cds.sam -as cds.fa -ar ref.sam -s nd.asm.fasta -n align1.anchors -R 1 -Q 1 -ns  -o 7a.maf
python3 closeGapUsingAnotherAssembly.py -r chr.7a.fa -q nd.asm.fasta -a align1.anchors -o chr.7a.filled.fa
```

Usage:

```
usage: closeGapUsingAnotherAssembly.py [-h] -r REF_FASTA_FILE -q QUERY_FASTA_FILE -a ANCHOR_FILE -o OUT_FILE

Close chromosome level assembly gap using the contig assembly from another assembly software or data

optional arguments:
  -h, --help            show this help message and exit
  -r REF_FASTA_FILE, --ref REF_FASTA_FILE
                        The reference genome of chromosome assembly in FASTA format
  -q QUERY_FASTA_FILE, --query QUERY_FASTA_FILE
                        The query genome of contig or scarfold or chromosome assembly in FASTA format
  -a ANCHOR_FILE, --anchor ANCHOR_FILE
                        The anchor file generated using AnchorWave
  -o OUT_FILE, --output OUT_FILE
                        The output File in FASTA format
```

Questions should go to songbaoxing168@163.com.
