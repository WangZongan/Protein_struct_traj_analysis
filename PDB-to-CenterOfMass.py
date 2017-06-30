#!/usr/bin/env python
'''
Given pdb file(s), calculate the RMSD values using Bio.PDB. 
'''
__author__  = 'Wang Zongan'
__version__ = '2016-10-11'

import os
import sys
import string
import numpy as np
import mdtraj as md


def main():
    traj = md.load_pdb(sys.argv[1])
    com = md.compute_center_of_mass(traj)

    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+'.CoMz.dat','w') as f:
        for e in com[:,2]*10:
            print >> f, '%.3f' % e

if __name__ == '__main__':
    main()

