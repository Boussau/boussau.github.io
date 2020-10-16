import sys
import os

inname = sys.argv[1]

nuc = ['A', 'C', 'G', 'T']
comp = [dict(), dict(), dict()]

for pos in range(3):
    for n in nuc:
        comp[pos][n] = 0

with open(inname, 'r') as infile:
    fields = infile.readline().rstrip('\n').split()
    ntaxa = int(fields[0])
    nsite = int(fields[1])
    ncodon = nsite // 3
    for line in infile:
        fields = line.rstrip('\n').split()
        if len(fields) == 2:
            name = fields[0]
            seq = fields[1]
            for i in range(ncodon):
                c = seq[3*i:3*(i+1)]
                for pos in range(3):
                    if c[pos] in nuc:
                        comp[pos][c[pos]] = comp[pos][c[pos]] + 1

for pos in [2]:
# for pos in range(3):
    tot = 0
    for n in nuc:
        tot = tot + comp[pos][n]
    print("pos : ", pos)
    for n in nuc:
        print(n, comp[pos][n] / tot)
    print()




