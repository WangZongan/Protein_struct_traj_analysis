#!/usr/bin/env python
'''
add backbone H to pdb file - orig Glen Hocky, mod AA

Copied from Aashish, modified by Wang Zongan.
Date: 2016-01-25
'''

import sys
from Bio.PDB import *

import time
import StringIO
import os, math
import numpy as np

NHLength=1 #Angstrom

def add_hn(prot):
    normalize = lambda vec: vec/np.sqrt(np.sum(vec*vec,axis=-1))

    atomlist=[]
    for residue in prot.get_residues():
        n_atom  = residue['N']
        ca_atom = residue['CA']
        c_atom  = residue['C']
        atomlist.append([n_atom,ca_atom,c_atom])

    for resnum,residue in enumerate(prot.get_residues()):
        if resnum==0 or residue.get_resname()=="PRO": continue

        c_atom  = atomlist[resnum-1][2]
        n_atom  = atomlist[resnum][0]
        ca_atom = atomlist[resnum][1]

        n_ca_vec = ca_atom.get_coord() - n_atom.get_coord()
        n_c_vec  = c_atom.get_coord()  - n_atom.get_coord()

        n_ca_vec = normalize(n_ca_vec)
        n_c_vec  = normalize(n_c_vec)

        sum_vec = n_ca_vec + n_c_vec
        sum_vec = normalize(sum_vec)

        NH_vec  = -sum_vec * NHLength + n_atom.get_coord()
        nh_atom = Atom.Atom('HN',NH_vec,1,1," "," HN ",1)
        residue.add(nh_atom)

    return prot


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Program used for adding missing HN atoms to pdb files.',
        usage ='use "python %(prog)s --help or -h" for more information')
    parser.add_argument('inputpdb', help='[required] input pdb file')
    args = parser.parse_args()

    protein = PDBParser(QUIET=True).get_structure('x',args.inputpdb)
    prot = add_hn(protein)

    io = PDBIO()
    io.set_structure(prot)
    io.save(os.path.splitext(os.path.basename(args.inputpdb))[0] + '-wH.pdb')


if __name__ == '__main__':
    main()

