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

__author__ = 'Wang.Zongan'
__version__ = '2016.05.13'

'''
Adapted from http://mdtraj.org/latest/examples/native-contact.html
'''
from itertools import combinations
def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`

    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')

    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])

    # compute the distances between these pairs in the native state
    def compute_distances(trajectory, pairs):
        '''
        return: distances, np.ndarray, shape=(n_frames, num_pairs), dtype=float
        '''
        return np.sqrt(np.sum(np.diff(trajectory.xyz[:,pairs],axis=2)**2,axis=-1))[:,:,0]

    heavy_pairs_distances = compute_distances(native[0], heavy_pairs)[0]

    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))

    # now compute these distances for the whole trajectory
    r = compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q


def calc_Qw(traj, native):
    """
    PNAS 111 (2014) 11031
    
    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used
        
    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of traj
    """
    
    nres = traj.n_residues

    ALPHA_CONST = 0.15
    A_CONST = 0.1 # nm
    
    # get the indices of all of the CA atoms
    ca = native.topology.select('name CA')
    ca_pairs = np.array([(i,j) for (i,j) in combinations(ca, 2) 
                         if abs(native.topology.atom(i).residue.index - 
                                native.topology.atom(j).residue.index) > 2])
    
    res_pairs = np.array([(i,j) for (i,j) in combinations(range(nres),2) if np.absolute(j-i) > 2 ])
    res_pairs_coefficient = (2*(A_CONST*np.absolute(np.diff(res_pairs, axis=1))**ALPHA_CONST))[:,0] # shape: (len(res_pair),)
        
    # compute the distances between these pairs in the native state
    def compute_distances(trajectory, pairs):
        '''
        return: distances, np.ndarray, shape=(n_frames, num_pairs), dtype=float
        '''
        return np.sqrt(np.sum(np.diff(trajectory.xyz[:,pairs],axis=2)**2,axis=-1))[:,:,0]
    
    # now compute these distances for the whole trajectory
    r = compute_distances(traj, ca_pairs)
    # and recompute them for just the native state
    r0 = compute_distances(native[0], ca_pairs)

    return np.sum(np.exp( -(r-r0)**2/res_pairs_coefficient**2), axis=1)*2./((nres-2)*(nres-3))


def calc_Qc(traj, native):
    """
    PNAS 111 (2014) 11031
    
    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used
        
    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of traj
    """
    
    nres = traj.n_residues

    ALPHA_CONST = 0.15
    A_CONST = 0.2 # nm
    NATIVE_CONTACT_CUTOFF = 0.95 # nm 
    
    # get the indices of all of the CA atoms
    ca = native.topology.select('name CA')
    ca_pairs = np.array([(i,j) for (i,j) in combinations(ca, 2) 
                         if abs(native.topology.atom(i).residue.index - 
                                native.topology.atom(j).residue.index) > 9])
    
    res_pairs = np.array([(i,j) for (i,j) in combinations(range(nres),2) if np.absolute(j-i) > 9 ])
        
    # compute the distances between these pairs in the native state
    def compute_distances(trajectory, pairs):
        '''
        return: distances, np.ndarray, shape=(n_frames, num_pairs), dtype=float
        '''
        return np.sqrt(np.sum(np.diff(trajectory.xyz[:,pairs],axis=2)**2,axis=-1))[:,:,0]
    
    ca_pairs_distances = compute_distances(native[0], ca_pairs)[0]
    native_contacts = ca_pairs[ca_pairs_distances < NATIVE_CONTACT_CUTOFF]
    print("Number of native contacts", len(native_contacts))
    
    # now compute these distances for the whole trajectory
    r = compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = compute_distances(native[0], native_contacts)
    
    res_pairs_coefficient = (
            2*(A_CONST*np.absolute(np.diff(res_pairs[ca_pairs_distances < NATIVE_CONTACT_CUTOFF], axis=1))**ALPHA_CONST)
        )[:,0] # shape: (len(res_pair),)

    return np.mean(np.exp( -(r-r0)**2 / res_pairs_coefficient**2), axis=1)


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
	return md.rmsd(traj, ref, ref_frame, atom_indices=atom_indices) 
    else:
        return md.rmsd(traj, ref, ref_frame)


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


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Calculate RMSD, contact score Q, and output potential.',
        usage ='use "%(prog)s --help" for more information')

    parser.add_argument('h5_file', help='[required] input H5 file')  
    parser.add_argument('pdb_traj_file', help='[required] input PDB trajectory file')
    parser.add_argument('pdb_ref_file', help='[required] input PDB references file')
    parser.add_argument('--stride', default=1, type=int,
        help = 'Output the potential by every number of stride.')  

    parser.add_argument('--rmsd', default=False, action='store_true', help = 'If turned on, calculate RMSD.') 
    parser.add_argument('--rmsd-selection', default=[], action='append', type=parse_segments, 
        help = 'If provided, only calculate RMSD for the selection.') 

    parser.add_argument('--q', default=False, action='store_true', help = 'If turned on, calculate Q.') 
    parser.add_argument('--qw', default=False, action='store_true', help = 'If turned on, calculate Qw.') 
    parser.add_argument('--qc', default=False, action='store_true', help = 'If turned on, calculate Qc.')
    parser.add_argument('--potential', default=False, action='store_true', help = 'If turned on, output the potential.')
    args   = parser.parse_args()

    h5file = args.h5_file
    trajfile = args.pdb_traj_file
    reffile = args.pdb_ref_file
    stride = args.stride 

    if args.rmsd or args.q or args.qw or args.qc:
        traj = md.load_pdb(trajfile)
        ref = md.load_pdb(reffile)

    if args.rmsd:
        rmsd = output_rmsd(traj, ref)
        with open(os.path.basename(h5file)+'.rmsd.dat','w') as f:
            for a in rmsd:
                f.write('%8.3f\n' % a)
        if args.rmsd_selection:
            rmsd = output_rmsd(traj, ref, selection=args.rmsd_selection)
            with open(os.path.basename(h5file)+'.rmsd_sel.dat','w') as f:
                for a in rmsd:
                    f.write('%8.3f\n' % a) 
    if args.q:
        q = best_hummer_q(traj, ref)
        with open(os.path.basename(h5file)+'.Q.dat','w') as f: 
            for a in q:
                f.write('%8.3f\n' % a) 
    if args.qw:
        qw = calc_Qw(traj, ref)
        with open(os.path.basename(h5file)+'.Qw.dat','w') as f:
            for a in qw: 
                f.write('%8.3f\n' % a) 
    if args.qc:
        qc = calc_Qc(traj, ref) 
        with open(os.path.basename(h5file)+'.Qc.dat','w') as f: 
            for a in qc: 
                f.write('%8.3f\n' % a) 
    if args.potential:
        pot = output_potential(h5file, stride)
        with open(os.path.basename(h5file)+'.potential.dat','w') as f:
            for a in pot:
                f.write('%.3f\n' % a) 

if __name__ == '__main__':
    main()


