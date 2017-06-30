#!/usr/bin/env python

import sys

contact_file = sys.argv[1]

with open(contact_file,'r') as f:
    lines = f.readlines()[1:]

r1 = []
r2 = []
for l in lines:
    l = l.strip()
    l = l.split()
    r1.append(l[0])
    r2.append(l[1])

restraint_string = ''
for i in range(len(r1)):
    restraint_string += '--restraint-group=%s,%s '%(r1[i],r2[i])

print restraint_string
