#!/usr/bin/env python
'''
Given pdb trajectory, calculate Rg.  
'''
__author__  = 'Wang Zongan'

import os
import sys
import string
import numpy as np
import cPickle as cp
import mdtraj as md



def main():
    traj = md.load(sys.argv[1])
    rg  = md.compute_rg(traj)*10  # unit: nm --> A

    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+'.Rg.dat', 'w') as f:
        for val in rg:
            print >> f, '%8.3f' % (val)


if __name__ == '__main__':
    main()
