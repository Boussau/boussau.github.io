import sys
import os

inname = sys.argv[1]
outname = sys.argv[2]

nuc = ['A', 'C', 'G', 'T']
with open(inname, 'r') as infile:
    with open(outname, 'w') as outfile:
        for i in range(7):
            line = infile.readline()
            outfile.write(line)
        for line in infile:
            fields = line.rstrip('\n').split()
            if len(fields) == 2:
                name = fields[0]
                seq = fields[1]
                nsite = len(seq)
                ncodon = nsite // 3
                seq2 = ['???' for i in range(ncodon)]
                for i in range(ncodon):
                    c = seq[3*i:3*(i+1)]
                    d = "???"
                    if c != "TAA" and c != "TAG" and c != "TGA" and c[0] in nuc and c[1] in nuc and c[2] in nuc:
                        d = c
                    seq2[i] = d
                outfile.write("{0}\t{1}\n".format(name, "".join(seq2)))
            else:
                outfile.write(line)




