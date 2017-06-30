# coding: utf-8
'''
'''
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np

def main():
    pklfile = sys.argv[1]
    flexible = int(sys.argv[2])

    with open(pklfile,'r') as f:
        bounds = cp.load(f)
    config = ''
    for bound in bounds:
        lb, ub = bound
        if ub - lb > 10:
            config += '--restraint-group=%i-%i '%(lb+5, ub-5)
        else:
            continue

    outputname = os.path.splitext(os.path.basename(pklfile))[0] + '.restraint_grp.'+str(flexible)+'flexible.config'
    with open(outputname,'w') as f:
        f.write(config)


if __name__ == '__main__':
    main()
