#!/usr/bin/env python
'''
This program is to obtain the positions of certain residues in the trajectory. 
'''

__author__  = 'Wang Zongan'
__version__ = '2017-03-28'

import os
import sys 
import string
import numpy as np
import cPickle as cp
from itertools import combinations

import Bio.PDB

import mdtraj.core.element as el
import mdtraj as md
from mdtraj.formats.registry import FormatRegistry
angstrom=0.1  # conversion to nanometer from angstrom

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *

##<-------------------- protein basics -------------------->##
backbone = ['C','CA','N','O']

base_sc_ref = {
    'ALA': np.array([-0.01648328,  1.50453228,  1.20193768]),
    'ARG': np.array([-0.27385093,  3.43874264,  2.24442499]),
    'ASN': np.array([-0.27119135,  2.28878532,  1.32214314]),
    'ASP': np.array([-0.19836569,  2.23864046,  1.36505725]),
    'CYS': np.array([-0.17532601,  1.92513503,  1.34296652]),
    'GLN': np.array([-0.28652696,  2.84800873,  1.60009894]),
    'GLU': np.array([-0.26377398,  2.80887008,  1.69621717]),
    'GLY': np.array([-1.56136239e-02, 5.46052464e-01, -5.67664281e-19]),
    'HIS': np.array([-0.32896151,  2.66635893,  1.42411271]),
    'ILE': np.array([-0.23956042,  2.26489309,  1.49776818]),
    'LEU': np.array([-0.23949426,  2.67123263,  1.3032201 ]),
    'LYS': np.array([-0.26626635,  3.18256448,  1.85836641]),
    'MET': np.array([-0.21000946,  2.79544428,  1.52568726]),
    'PHE': np.array([-0.27214755,  2.83761534,  1.45094383]),
    'PRO': np.array([-1.10993493,  0.89959734,  1.41005877]),
    'SER': np.array([-0.00692474,  1.56683138,  1.475341  ]),
    'THR': np.array([-0.14662723,  1.80061252,  1.42785569]),
    'TRP': np.array([-0.01433503,  3.07506159,  1.56167948]),
    'TYR': np.array([-0.2841611 ,  3.02555746,  1.50123341]),
    'VAL': np.array([-0.02436993,  1.97251406,  1.32782961])}

model_geom = np.zeros((3,3))
model_geom[0] = (-1.19280531, -0.83127186, 0.)  # N
model_geom[1] = ( 0.,          0.,         0.)  # CA
model_geom[2] = ( 1.25222632, -0.87268266, 0.)  # C
model_geom -= model_geom.mean(axis=0)

three_letter_aa = dict(
    A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE', K='LYS', L='LEU', 
    M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER', T='THR', V='VAL', W='TRP', Y='TYR')

aa_num = dict([(k,i) for i,k in enumerate(sorted(three_letter_aa.values()))])
one_letter_aa = dict([(v,k) for k,v in three_letter_aa.items()])


def rmsd_transform(target, model):
    assert target.shape == model.shape == (model.shape[0], 3)
    base_shift_target = target.mean(axis=0)
    base_shift_model  = model .mean(axis=0)

    target = target - target.mean(axis=0)
    model  = model  - model .mean(axis=0)

    R = np.dot(target.T, model)
    U,S,Vt = np.linalg.svd(R)
    if np.linalg.det(np.dot(U,Vt))<0.:
        Vt[:,-1] *= -1.  # fix improper rotation
    rot   = np.dot(U,Vt)
    shift = base_shift_target - np.dot(rot, base_shift_model)
    return rot, shift


def calculate_sidechain_CM(N, CA, C, resType):
    '''
    N, CA, C : coordinates of N, CA, C atoms of one given residue
    resType  : 3-letter aa type
    '''
    assert N.shape == CA.shape == C.shape == (3, )
    z = np.vstack((N, CA, C))
    #z = np.concatenate((N, CA, C), axis=0)
    rot, trans = rmsd_transform(z, model_geom)
    return np.dot(base_sc_ref[resType],rot.T) + trans


##<-------------------- John Jumper's functions for loading h5 file -------------------->##
## adpated so that I don't need to import mdtraj
def _output_groups(t):
    i=0
    while 'output_previous_%i'%i in t.root:
        yield t.get_node('/output_previous_%i'%i)
        i += 1
    if 'output' in t.root:
        yield t.get_node('/output')
        i += 1


def load_upside_traj(fname, sel='CM', stride=1, target_pos_only=False):
    import tables as tb
    with tb.open_file(fname) as t:
        start_frame = 0
        total_frames_produced = 0
        xyz = []
        if target_pos_only:
            xyz.append(t.root.target.pos[:,:,0])
            total_frames_produced = 1
            start_frame=1
        else:
            for g_no, g in enumerate(_output_groups(t)):
                # take into account that the first frame of each pos is the same as the last frame before restart
                # attempt to land on the stride
                sl = slice(start_frame,None,stride)
                xyz.append(g.pos[sl,0])
                total_frames_produced += g.pos.shape[0]-(1 if g_no else 0)  # correct for first frame
                start_frame = 1 + stride*(total_frames_produced%stride>0) - total_frames_produced%stride
        xyz = np.concatenate(xyz,axis=0)  # N, CA, C, (n_frame, 3*n_res, 3)
        seq = t.root.input.sequence[:]
    print "Finished loading upside trajectory."
    
    nmodel, nres = xyz.shape[0], len(seq)
    assert xyz.shape == (nmodel, nres*3, 3)
    coords = np.zeros((nmodel, nres, 3))
    if sel == 'CA':
        for nm in range(nmodel):
            for nr in range(nres):
                coords[nm][nr] = xyz[nm][nr*3+1]
        return coords.astype('f4') 
    elif sel == 'CM':
        print "Now, use the static reference sidechain CM positions implied by N, CA, and C."
        for nm in range(nmodel):
            for nr in range(nres):
                coords[nm][nr] = calculate_sidechain_CM(xyz[nm][nr*3], xyz[nm][nr*3+1], xyz[nm][nr*3+2], seq[nr])
        print "Finished computing all referred sidechain positions."
        return coords.astype('f4')
    else:
        raise ValueError('--contact-type must be either CM or CA')


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


def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Obtain the postions of certain residues from the h5 file. ',
        usage ='use "%(prog)s --help" for more information')
    parser.add_argument('h5', help='[required] input Upside-output H5 file')
    parser.add_argument('atomtype', default='CA', type=str, help='[required] atom type')
 
    parser.add_argument('--residue-group', default=[], action='append', type=parse_segments,
        help = 'Multiple residue groups may be specified by giving the --residue-group flag multiple times. ' +
               'Note: residue indices start from 0.' )

    parser.add_argument('--output-fname', default=None, type=str, help='Output file name.')
    args = parser.parse_args()

    # atoms selection & obtain coordinates
    coords       = load_upside_traj(args.h5, args.atomtype)
    nmodel, nres = coords.shape[0], coords.shape[1]
    print nmodel, nres

    if args.output_fname is not None:
        for i, rg in enumerate(args.residue_group):
            with open(args.output_fname+'.%i.pkl'%i, 'w') as f:
                cp.dump((rg, coords[:, rg]), f, -1)
                print coords[:, rg].shape
        
    return True


if __name__ == '__main__':
    from time import time
    sta = time()
    main()
    print 'running time: %s' % sec_to_hr_min_sec(time() - sta)




