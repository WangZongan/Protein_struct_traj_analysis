#!/usr/bin/env python
'''
Given pdb file(s), calculate the RMSD matrix. 
'''
__author__  = 'Wang Zongan'
__version__ = '2016-08-05'

import os
import sys
import string
import numpy as np
import cPickle as cp
import mdtraj as md

def calculate_rmsd_mat_mdtraj(com_traj, ref_traj):
    from msmbuilder.featurizer import RMSDFeaturizer
    return RMSDFeaturizer(ref_traj).partial_transform(com_traj)

def calculate_drid_mat(traj):
    from msmbuilder.featurizer import DRIDFeaturizer
    return DRIDFeaturizer().partial_transform(traj)

def calculate_dihe_mat(traj):
    from msmbuilder.featurizer import DihedralFeaturizer
    return DihedralFeaturizer().partial_transform(traj)

def calculate_contact_mat(traj,scheme):
    from msmbuilder.featurizer import ContactFeaturizer
    return ContactFeaturizer(scheme=scheme).partial_transform(traj)

def calculate_tmscore_mat(traj):
    pass 


def main():
    import argparse
    parser = argparse.ArgumentParser(
        usage ='use "python %(prog)s --help" for more information',
        description = 'Featurize the given pdb trajectory into a vectorizable space. ')
    parser.add_argument('pdb', help='[required] input PDB file.')
    parser.add_argument('--reference-pdb', type=str, 
        help = "If the reference pdb is not provided, " +  
               "the first model in the pdb structure will be used as the reference, if necessary. ")
    parser.add_argument('--featurizer', type=str, default='RMSD',
        help = "Choices of featurizer: 'RMSD', 'DRID', 'dihedral', 'contact'.")
    parser.add_argument('--output', type=str, default='output', 
        help = 'Output the pkl file of the shape (n_samples, n_features).')
    args = parser.parse_args()

    traj = md.load(args.pdb)
    if args.featurizer == 'RMSD': 
        if args.reference_pdb is not None:
            ref_traj = md.load(args.reference_pdb)
            rmsd_mat = calculate_rmsd_mat_mdtraj(traj, ref_traj)
        else:
            'Note: no reference trajectory provided, use the input trajectory itself as reference.'
            rmsd_mat = calculate_rmsd_mat_mdtraj(traj, traj)
        with open(args.output+'.RMSD.pkl','w') as f:
            cp.dump(rmsd_mat, f, -1)
    elif args.featurizer == 'DRID': 
        with open(args.output+'.DRID.pkl','w') as f:
            cp.dump(calculate_drid_mat(traj), f, -1)
    elif args.featurizer == 'dihedral':
        with open(args.output+'.dihedral.pkl','w') as f:
            cp.dump(calculate_dihe_mat(traj), f, -1)
    elif args.featurizer == 'contact':
        with open(args.output+'.CA_contact.pkl','w') as f: 
            cp.dump(calculate_ca_contact_mat(traj), f -1)


if __name__ == '__main__':
    main()

