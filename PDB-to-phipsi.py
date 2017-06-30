#!/usr/bin/env python
'''
Generate a text file that contains information for a contact energy function for UPSIDE.

The first line of the file should be a header: "index angle_type start end width energy".
The remaining lines should contain space separated values.

good_res_num = (f[1]=='phi' and 0 < res_num <= n_res-1) or (f[1]=='psi' and 0 <= res_num < n_res-1) 

Example file:
index angle_type start end width energy
0 psi 
1 phi
1 psi 
2
....

start, end, width are all in the unit of degree. 

The formula of the contact energy:

                                   energy                            1
    dihedral_energy = ------------------------------- * --------------------------- ,
                                 1                                 1
                      1 + exp[ ----- * (- x + start)]   1 + exp[ ----- * (x - end)]
                               width                             width

    in which x is the dihedral angle.
'''

__author__  = 'Wang Zongan'
__version__ = '2016-06-22'

import string
import numpy as np
import argparse
import mdtraj as md
import Bio.PDB 

three_letter_aa = dict(
        A='ALA', C='CYS', D='ASP', E='GLU',
        F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN',
        P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

aa_num = dict([(k,i) for i,k in enumerate(sorted(three_letter_aa.values()))])

one_letter_aa = dict([(v,k) for k,v in three_letter_aa.items()])

def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
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

def calc_phi(traj):
    # phi: 1 ~ nres-1
    _, phi = md.compute_phi(traj) # phi in rand
    return phi[0]*180/np.pi

def calc_psi(traj):
    # psi: 0 ~ nres-2
    _, psi = md.compute_psi(traj)
    return psi[0]*180/np.pi

def calc_ss(traj):
    return md.compute_dssp(traj) 

larger  = lambda a,b: a if a>=b else b
smaller = lambda a,b: a if a<=b else b

def main():
    parser = argparse.ArgumentParser(description='Prepare input file',
        usage ="use 'python %(prog)s --help' for more information")
    parser.add_argument('pdb',    help='[required] input .pdb file')
    parser.add_argument('output', help='[required] output file')
    parser.add_argument('--width' , default=5. , type=float,
        help = 'Parameter of width. Defaulf value is 5.')
    parser.add_argument('--energy', default=-3., type=float,
        help = 'Parameter of energy. Defaulf value is -3.')
    parser.add_argument('--flexible', default=10. , type=float,
        help = 'Parameter for start = x - flexible, end = x + flexible. Defaulf value is 10.')
    parser.add_argument('--restraint-group', default=[], action='append', type=parse_segments,
        help = 'If provided, only the residues in the restaint groups will be included in the contact energy calculation. ' +
               'Multiple residue groups may be specified by giving the --restraint-group flag multiple times ' +
               'with different filenames. Note: residue indices starts from 0.' )
    args = parser.parse_args()

    # Note: only calculate the dihedral angles for the first model in the pdb
    model   = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', args.pdb)[0] 
    nchains = len([chain for chain in model])
    nres    = len([res for res in model.get_residues() if res.get_id()[0] == ' '])

    dihedrals = []
    for chain in model :
        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
        dihedrals.extend(poly.get_phi_psi_list())
    phi = np.array([d[0]*180/np.pi if d[0] != None else None for d in dihedrals]) 
    psi = np.array([d[1]*180/np.pi if d[1] != None else None for d in dihedrals]) 

    with open(args.output,'w') as f:
        print >> f, 'index angle_type start end width energy'
        if args.restraint_group:
            indices = []
            for rg in args.restraint_group:
                indices.extend(rg)
            print 'Indices selected: %s' % indices
        else:
            print 'All residue are selected for calculation.'
            indices = np.arange(nres)

        for idx in indices:
            if phi[idx] != None:
                print >> f, '%i %s %.3f %.3f %.3f %.3f' % (
                        idx,
                        'phi',
                        larger (phi[idx]-args.flexible,-179.999),
                        smaller(phi[idx]+args.flexible, 179.999),
                        args.width,
                        args.energy)
            if psi[idx] != None:
                print >> f, '%i %s %.3f %.3f %.3f %.3f' % (
                        idx,
                        'psi',
                        larger (psi[idx]-args.flexible,-179.999),
                        smaller(psi[idx]+args.flexible, 179.999),
                        args.width,
                        args.energy)


if __name__ == '__main__':
    main()
