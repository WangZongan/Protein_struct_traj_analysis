#!/usr/bin/env python
import os
import sys
import string
import numpy as np
import mdtraj as md


def calculate_ss(pdbfilename):
    prot = md.load(pdbfilename)
    ss   = md.compute_dssp(prot)
    secseq = ''.join((ele for ele in ss[0]))
    return secseq


def get_basename_no_ext(path):
    return os.path.splitext(os.path.basename(path))[0]

def main():
    native = sys.argv[1]
    compare = sys.argv[2]

    ss_n = calculate_ss(native)
    ss_c = calculate_ss(compare)

    assert len(ss_n) == len(ss_c)
    with open('%s.%s.compare_ss.dat'%(get_basename_no_ext(native),get_basename_no_ext(compare)),'w') as f:
        print >> f, ss_n
        print >> f, ss_c
        
        percent_same = 0
        for i in range(len(ss_n)):
            if ss_n[i] == ss_c[i]:
                percent_same += 1
        print >> f, '%.3f' % (percent_same*1. / len(ss_n)*1.)


if __name__ == '__main__':
    main()



