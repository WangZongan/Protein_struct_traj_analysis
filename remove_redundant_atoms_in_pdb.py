#!/usr/bin/env python
'''
Remove redundant atoms in PDB files. 

Wang Zongan.
Date: 2016-11-28
'''
import string, math, sys, os
import numpy as np

'''
0         1         2         3         4         5         6         7     
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM    208  OE1AGLU A  31      26.496  11.602  18.389  0.50 29.30           O  
ATOM    209  OE1BGLU A  31      28.714  11.470  17.870  0.50 23.19           O 

ATOM      1  N   GLY A  -1      32.647   1.960 -18.519  1.00 69.77      A    N
'''

def main():
    pdbfile = sys.argv[1]
    
    lines = []
    with open(pdbfile,'r') as f:
        for l in f.readlines():
            if l[0:4] == 'ATOM':
                if int(l[22:26].strip()) >= 0:
                    lines.append(l)
    print len(lines)

    i = 0 
    j = 1 
    bad_lines_idx = []
    while i < len(lines)-1:
        if lines[i][13:16] == lines[i+j][13:16] and lines[i][17:20] == lines[i+j][17:20] and lines[i][21] == lines[i+j][21]:
            bad_lines_idx.append(i+j)
            j += 1
        else:
            i = i+j
            j = 1

    with open(os.path.splitext(os.path.basename(pdbfile))[0]+'protein.pdb', 'w') as f:
        for il, l in enumerate(lines):
            if il not in bad_lines_idx:
                if l[16] != ' ':
                    new_l = l[0:16]+' '+l[17:]
                    print >> f, new_l,
                else:
                    print >> f, l, 

if __name__ == '__main__':
    main()

