#!/usr/bin/env python
'''
Convert vtf format trajectories to pdb format (or pdt format).
'''

__author__  = 'Wang Zongan'
__version__ = '2016-10-27'

import os
import string
import numpy as np

from collections import namedtuple
atom = namedtuple('atom','atomid atom_name resid res_name chainid trajectory'.split())


def read_vtf(vtf_file, CAonly=False, CPRasPRO=False):
    '''
    VTF
    -------------------------------------------------
    atom 0 name N  resid 0 resname GLY segid s0
    atom 1 name CA resid 0 resname GLY segid s0
    atom 2 name C  resid 0 resname GLY segid s0
    bond 0:1
    bond 1:2
    atom 3 name O  resid 0 resname GLY segid s0
    bond 2:3
    ...

    timestep ordered
    -13.447 2.052 28.000
    -12.337 1.447 29.126
    -11.131 0.555 28.323
    -10.080 0.221 28.889
    ...
    '''
    atoms = []

    # determine number of frames
    os.system("grep -n timestep %s > timestep.temp" % vtf_file)
    with open('timestep.temp','r') as f:
        lines = f.readlines()
        nframe = len(lines)
        timestep_idx = []
        for line in lines:
            line = line.strip().split()[0].split(':')
            timestep_idx.append(int(line[0]))
    print 'nframe: %i' % nframe

    # read atoms
    os.system("grep atom %s > atom.temp" % vtf_file)
    with open('atom.temp','r') as f:
        lines = f.readlines()
        natom = len(lines)
        for line in lines:
            line = line.strip().split()
            atomid, atomname, resid, resname = int(line[1]), line[3], int(line[5]), line[7]
            atom_ = atom(atomid, atomname, resid, resname, 'A', np.zeros((nframe,3)))
            atoms.append(atom_)
    print 'natom : %i' % len(atoms)
    
    # delete temp files
    os.system("rm timestep.temp atom.temp")

    # read coordinates
    with open(vtf_file,'r') as f:
        lines = f.readlines()
        for i in range(nframe):
            sublines = lines[timestep_idx[i]:timestep_idx[i]+natom]
            for j,atom_ in enumerate(atoms):
                line = sublines[j].strip().split()
                atom_.trajectory[i] = float(line[0]),float(line[1]),float(line[2])
    
    if CAonly:
        CA = [atom_ for atom_ in atoms if atom_.atom_name == 'CA']
        atom_count = 0
        for atom_ in CA:
            atom_ = atom_._replace(atomid = atom_count)
            atom_count += 1
        return CA
    else:
        # locate CPR atoms 
        '''
        Example: 
        atom 313 name N  resid 63 resname CPR segid s0
        atom 314 name CA resid 63 resname CPR segid s0
        atom 315 name C  resid 63 resname CPR segid s0
        bond 310:313
        bond 313:314
        bond 314:315
        atom 316 name H  resid 63 resname CPR segid s0
        atom 317 name O  resid 63 resname CPR segid s0

        Change CRP to PRO; skip H. 
        '''
        if CPRasPRO:
            CPR_H = [atom_ for atom_ in atoms if atom_.res_name == 'CPR' and atom_.atom_name == 'H']
            nCPR_H = len(CPR_H)
            atom_count = 0
            new_atoms  = []
            for atom_ in atoms:
                if atom_.res_name != 'CPR':
                    atom_ = atom_._replace(atomid = atom_count)
                    new_atoms.append(atom_)
                    atom_count += 1
                elif atom_.res_name == 'CPR' and atom_.atom_name != 'H':
                    atom_ = atom_._replace(atomid  = atom_count)
                    atom_ = atom_._replace(res_name= 'PRO')
                    new_atoms.append(atom_)
                    atom_count += 1
            assert (nCPR_H + len(new_atoms)) == natom
            return new_atoms
        else:
            return atoms


def write_pdb(fname, atoms, chain_break):
    f = open(fname,'w')

    '''
    PDB
    --------------------------------------------------------------------------------
    0         1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ATOM      5  N   LEU X   1       8.802  -8.387  31.340  0.00  0.00      s0
    ATOM      6  CA  LEU X   1       8.643  -6.895  31.816  0.00  0.00      s0
    ATOM      7  C   LEU X   1       8.622  -5.988  30.632  0.00  0.00      s0
    ATOM      8  H   LEU X   1       9.565  -8.783  31.528  0.00  0.00      s0
    ATOM      9  O   LEU X   1       7.787  -5.072  30.662  0.00  0.00      s0
    ...
    '''

    f.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n')
    natom = len(atoms)
    print 'N saved atoms / frame : %i' % natom
    nframe = atoms[0].trajectory.shape[0]

    resid_max = max([atom_.resid for atom_ in atoms])
    new_chain_break = [0]
    if chain_break  == None:
        nchain = 1
    else:
        nchain = len(chain_break)+1
        new_chain_break.extend(chain_break)
    new_chain_break.append(resid_max+1)
    chain_nm = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[0:nchain]

    atom_ch_id = []
    for atom_ in atoms:
        for i in range(nchain):
            if atom_.resid >= new_chain_break[i] and atom_.resid < new_chain_break[i+1]:
                atom_ch_id.append(chain_nm[i])

    for i in range(nframe):
        for j,atom_ in enumerate(atoms):
            f.write('ATOM  %5i  %-3s%4s %s%4i    %8.3f%8.3f%8.3f  0.00  0.00      s0\n' %
                (atom_.atomid,
                 atom_.atom_name,
                 atom_.res_name,
                 atom_ch_id[j],
                 atom_.resid,
                 atom_.trajectory[i][0],
                 atom_.trajectory[i][1],
                 atom_.trajectory[i][2]))
        f.write('ENDMDL\n')


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Prepare input file',
        usage ='use "%(prog)s --help" for more information')
    parser.add_argument('vtf', help='[required] input VTF file')
    parser.add_argument('pdb', help='[required] output PDB file')
    parser.add_argument('--chain-break', default=None,
        help = 'if there are chain breaks in the protein or suppose the structure has multiple chains, '+
               'with the number(s) of residue(s) provided, the chains will denoted as A, B, etc.' +
               'Example form: 20,50,100 (no spacing allowed). ' + 
               'Chain breaks are the first residues of each chain; 0 need not be provided.')
    parser.add_argument('--stride', default=1, type=int,
        help = 'save the trajectory by every number of stride.')
    parser.add_argument('--CPRasPRO', default=False, action='store_true',
        help = 'If turned on, not write CPR, instead write PRO still.')
    parser.add_argument('--CAonly', default=False, action='store_true',
        help = 'If turned on, only write out CA atoms. Helpful when you want to calculate CA-RMSD.')
    args = parser.parse_args()

    chain_break = args.chain_break
    if chain_break != None:
        chain_break = [int(cb) for cb in chain_break.strip().split(',')]

    stride = args.stride
    atoms = read_vtf(args.vtf, CAonly=args.CAonly, CPRasPRO=args.CPRasPRO)

    if stride == 1:
        write_pdb(args.pdb, atoms, chain_break)
    else:
        new_atoms = []
        for atom_ in atoms:
            new_traj = atom_.trajectory[::stride].astype('f4')
            new_atom_ = atom(atom_.atomid, atom_.atom_name, atom_.resid, atom_.res_name, atom_.chainid, new_traj)
            new_atoms.append(new_atom_)
        write_pdb(args.pdb, new_atoms, chain_break)

if __name__ == '__main__':
    main()
