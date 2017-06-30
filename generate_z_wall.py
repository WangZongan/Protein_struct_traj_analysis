#!/usr/bin/env python
import sys

thick = float(sys.argv[1])
nres = int(sys.argv[2])
output = sys.argv[3]

with open(output,'w') as f:
    print >> f, 'residue z0 radius spring_constant'
    for i in range(nres):
        print >> f, '%i 0 %.3f 100' % (i, thick)

