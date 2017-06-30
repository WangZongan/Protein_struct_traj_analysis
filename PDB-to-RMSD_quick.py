#!/usr/bin/env python
'''
Given pdb trajectory, perform FAST analysis.. 
'''
__author__  = 'Wang Zongan'
__version__ = '2016-10-27'

import os
import sys
import string
import numpy as np
import cPickle as cp
import mdtraj as md

def calculate_rmsd_mat_mdtraj(com_traj, ref_traj):
    from msmbuilder.featurizer import RMSDFeaturizer
    return RMSDFeaturizer(ref_traj).partial_transform(com_traj)


def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def main():
    traj = md.load(sys.argv[1])
    target = md.load(sys.argv[2])
    target_rmsd = calculate_rmsd_mat_mdtraj(traj, target)

    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+'.targetRMSD.dat', 'w') as f:
        for val in target_rmsd:
            print >> f, '%8.3f' % (val)


if __name__ == '__main__':
    from time import time
    sta = time()
    main()
    print 'running time: %s' % sec_to_hr_min_sec(time() - sta)
