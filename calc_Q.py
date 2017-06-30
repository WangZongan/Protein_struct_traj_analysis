#!/usr/bin/env python
import re,gc
import sys,os,time,math,string
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np
import tables as tb
import mdtraj as md

from itertools import combinations


__author__ = 'Wang.Zongan'
__version__ = '2017.05.16'


def calc_Q(traj, native, pairs):
    """
    PNAS 113 (2016) 2098
    
    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for.
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used.
    pairs: np.array, shape = (npairs, 2)
        Pairs of residue indices.
    
    Formula
    -------
          1                     -(r_ij-r0_ij)^2 
    Q = ----- * sum_(i,j) exp( ----------------- ) 
          Np                     2*sigma_ij^2 
    
    sigma_ij = |i-j|^0.15

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of traj.
    """
    
    nres   = traj.n_residues
    npairs = len(pairs)

    # get pair coefficients
    PAIR_COEFF = 0.15 
    pair_coefficient = 2*np.array([np.power(np.absolute(pairs[i][0] - pairs[i][1]), PAIR_COEFF)**2 for i in range(npairs)])

    # get the indices of all of the CA atoms
    ca = native.topology.select('name CA')
        
    # compute the distances between these pairs in the native state
    def compute_distances(trajectory, pairs):
        '''
        return: distances, np.ndarray, shape=(n_frames, num_pairs), dtype=float
        '''
        return np.sqrt(np.sum(np.diff(trajectory.xyz[:,pairs],axis=2)**2,axis=-1))[:,:,0]*10 # unit: nm --> A
    
    # now compute these distances for the whole trajectory
    r = compute_distances(traj, pairs)
    # and recompute them for just the native state
    r0 = compute_distances(native[0], pairs)

    return np.sum(np.exp(-(r-r0)**2/pair_coefficient), axis=1)/npairs


def output_potential(h5file, stride):
    config = tb.open_file(h5file)
    pot = config.root.output.potential[::stride]
    config.close()
    return pot


def output_rmsd(traj, ref, ref_frame=0, selection=None): 
    if selection is not None:
        print selection
        atom_indices = []
        for i in selection[0]:
            atom_indices.extend(traj.topology.select('residue %i and name CA' % i))
        atom_indices = np.array(atom_indices)
        return md.rmsd(traj, ref, ref_frame, atom_indices=atom_indices)*10  # unit: nm --> A
    else:
        return md.rmsd(traj, ref, ref_frame)*10  # unit: nm --> A


def get_res_pairs(res_indices, short_or_long_range='all'):
    '''
    res_indices: np.array, shape = (n, )
    short_or_long_range:
        Available choices: 
            all: |i-j| >= 3
            short: 3 <= |i-j| <= 8
            long: |i-j| >= 9
    '''
    if short_or_long_range == 'all':
        return np.array([(i,j) for (i,j) in combinations(res_indices, 2) if abs(i-j) >= 3])
    elif short_or_long_range == 'short':
        return np.array([(i,j) for (i,j) in combinations(res_indices, 2) if 3 <= abs(i-j) <= 8])
    elif short_or_long_range == 'long':
        return np.array([(i,j) for (i,j) in combinations(res_indices, 2) if abs(i-j) >= 9])
    else:
        raise ValueError('short_or_long_range should be all, short, or long.')


def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
    import argparse 
    import re

    if re.match('^([0-9]+(-[0-9]+)?)(,[0-9]+(-[0-9]+)?)*$', s) is None:
        raise argparse.ArgumentTypeError('segments must be of the form 10-30,45,72-76 or similar')

    def parse_seg(x):
        atoms = x.split('-')
        if len(atoms) == 1:
            return np.array([int(atoms[0])])
        elif len(atoms) == 2:
            return np.arange(int(atoms[0]),1+int(atoms[1]))  # inclusive on both ends
        else:
            raise RuntimeError('the impossible happened.  oops.')

    ints = np.concatenate([parse_seg(a) for a in s.split(',')])
    ints = np.array(sorted(set(ints)))   # remove duplicates and sort
    return ints


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0] 


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Calculate RMSD, contact score Q, and output potential.',
        usage ='use "%(prog)s --help" for more information')

    parser.add_argument('h5_file', help='[required] input H5 file')  
    parser.add_argument('pdb_traj_file', help='[required] input PDB trajectory file')
    parser.add_argument('pdb_ref_file', help='[required] input PDB references file')

    parser.add_argument('--rmsd', default=False, action='store_true', help = 'If turned on, calculate RMSD.') 
    parser.add_argument('--rmsd-selection', default=[], action='append', type=parse_segments, 
        help = 'If provided, only calculate RMSD for the selection.') 

    parser.add_argument('--q', default=False, action='store_true', help = 'If turned on, calculate Q.') 
    parser.add_argument('--residue-group', default=[], action='append', type=parse_segments, 
        help = 'residue group for residue pairs in contact.')
    parser.add_argument('--pair-range', type=str, default='all', help = 'Unique residue pairs are included.') 
    parser.add_argument('--output-midfix', type=str, default='all', help = 'midfix in the name of the output file of Q values.')

    parser.add_argument('--potential', default=False, action='store_true', help = 'If turned on, output the potential.')
    parser.add_argument('--stride', default=1, type=int, help = 'Output the potential by every number of stride.')  
    args   = parser.parse_args()

    h5file   = args.h5_file
    trajfile = args.pdb_traj_file
    reffile  = args.pdb_ref_file
    stride   = args.stride 

    if args.rmsd or args.q:
        traj = md.load_pdb(trajfile)
        ref  = md.load_pdb(reffile)
        nres = traj.n_residues

    if args.rmsd:
        rmsd = output_rmsd(traj, ref)
        with open(bsnm(h5file)+'.rmsd.dat','w') as f:
            for a in rmsd:
                f.write('%8.3f\n' % a)
        if args.rmsd_selection:
            rmsd = output_rmsd(traj, ref, selection=args.rmsd_selection)
            with open(bsnm(h5file)+'.rmsd_sel.dat','w') as f:
                for a in rmsd:
                    f.write('%8.3f\n' % a) 

    if args.q:
        if args.residue_group:
            res_ids = []
            for rg in args.residue_group:
                res_ids = list(set(res_ids).union(rg))
            res_ids = list(set(res_ids).intersection(np.arange(nres)))
            pairs = get_res_pairs(np.array(res_ids), args.pair_range)
        else:
            pairs = get_res_pairs(np.arange(nres), args.pair_range)
        print '%i residue pairs' % len(pairs)
        
        q = calc_Q(traj, ref, pairs)
        with open('%s.%s.Q.dat' % (bsnm(h5file), args.output_midfix) ,'w') as f: 
            for a in q:
                f.write('%8.3f\n' % a) 

    if args.potential:
        pot = output_potential(h5file, stride)
        with open(bsnm(h5file)+'.potential.dat','w') as f:
            for a in pot:
                f.write('%.3f\n' % a) 

if __name__ == '__main__':
    main()







