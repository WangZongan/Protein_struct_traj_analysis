#!/usr/bin/env python
import sys

n_systems = int(sys.argv[1])
assert n_systems > 1

s = ['%i-%i'%(i,i+1) for i in range(n_systems-1)]

with open('swp1','w') as f:
    print >> f, ','.join(s[0::2])
with open('swp2','w') as f:
    print >> f, ','.join(s[1::2])
