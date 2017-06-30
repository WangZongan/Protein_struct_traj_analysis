#!/usr/bin/env python
import os
import sys

with open(sys.argv[1],'r') as f:
    lines = f.readlines()
    fasta = lines[1].strip()

print sys.argv[1]
print len(fasta)
print fasta

half1 = fasta[:len(fasta)/2]
half2 = fasta[len(fasta)/2:]
new_fasta = half2 + half1

print half1 + half2
print new_fasta
assert len(new_fasta) == len(fasta)

with open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+'.invert_half.fasta','w') as f:
    print >> f, '>'
    print >> f, new_fasta
