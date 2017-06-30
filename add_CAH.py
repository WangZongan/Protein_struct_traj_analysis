#!/usr/bin/env python
'''
Add HN, CB, and Calpha H to pdb file.
For GLY, add both Calpha H's; for non-GLY's without Calpha H, add Calpha H. 

The original program for adding HN to pdb is given by Aashish. 
'''

__author__ = 'Wang Zongan'
__version__ = '2016.03.15'

import sys
import Bio.PDB

import time
import StringIO
import os,math
import numpy as np

# top_all27_prot_lipid.inp
# HSD for HIS 
N_NH_Length = { # Angstrom
    'ALA':0.9996, 'CYS':0.9982, 'ASP':0.9966, 'GLU':0.9961, 'PHE':0.9987, 
    'GLY':0.9992, 'HIS':0.9988, 'ILE':0.9978, 'LYS':0.9988, 'LEU':0.9979, 
    'MET':0.9978, 'ASN':0.9992, 'PRO':np.nan, 'GLN':0.9984, 'ARG':0.9973,  
    'SER':0.9999, 'THR':0.9995, 'VAL':0.9966, 'TRP':0.9972, 'TYR':0.9986 }

CA_CB_Length = { # Angstrom
    'ALA':1.5461, 'CYS':1.5584, 'ASP':1.5619, 'GLU':1.5516, 'PHE':1.5594, 
    'GLY':1.0817, 'HIS':1.5519, 'ILE':1.5681, 'LYS':1.5568, 'LEU':1.5543, 
    'MET':1.5546, 'ASN':1.5627, 'PRO':1.5399, 'GLN':1.5538, 'ARG':1.5552, 
    'SER':1.5585, 'THR':1.5693, 'VAL':1.5660, 'TRP':1.5560, 'TYR':1.5606 }

CA_H_Length = { # Angstrom
    'ALA':1.0840, 'CYS':1.0837, 'ASP':1.0841, 'GLU':1.0828, 'PHE':1.0832, 
    'GLY':1.0814, 'HIS':1.0830, 'ILE':1.0826, 'LYS':1.0833, 'LEU':1.0824, 
    'MET':1.0832, 'ASN':1.0848, 'PRO':1.0837, 'GLN':1.0832, 'ARG':1.0836, 
    'SER':1.0821, 'THR':1.0817, 'VAL':1.0828, 'TRP':1.0835, 'TYR':1.0833 }

normalize = lambda vec: vec/np.sqrt(np.sum(vec*vec,axis=-1))

def add_hn(prot):
    atomlist=[]
    for residue in prot.get_residues():
        n_atom  = residue['N']
        ca_atom = residue['CA']
        c_atom  = residue['C']
        atomlist.append([n_atom,ca_atom,c_atom])

    for resnum,residue in enumerate(prot.get_residues()):
        if resnum==0 or residue.get_resname()=="PRO": continue
        if 'HN' in [atom.get_name() for atom in residue.get_atom()]: continue

        c_atom  = atomlist[resnum-1][2]
        n_atom  = atomlist[resnum][0]
        ca_atom = atomlist[resnum][1]

        n_ca_vec = ca_atom.get_coord() - n_atom.get_coord()
        n_c_vec  = c_atom.get_coord()  - n_atom.get_coord()

        n_ca_vec = normalize(n_ca_vec)
        n_c_vec  = normalize(n_c_vec)

        sum_vec = n_ca_vec + n_c_vec
        sum_vec = normalize(sum_vec)

        NH_vec  = -sum_vec * N_NH_Length[residue.get_resname()] + n_atom.get_coord()
        nh_atom = Bio.PDB.Atom.Atom('HN',NH_vec,1,1," "," HN ",1)
        residue.add(nh_atom)

    return prot

def add_CB(prot):
    for resnum, residue in enumerate(prot.get_residues()):
        # get atom coordinates as vectors
        n_vec  = residue['N'].get_vector()
        c_vec  = residue['C'].get_vector() 
        ca_vec = residue['CA'].get_vector()
        # center at origin
        n_vec = n_vec - ca_vec
        c_vec = c_vec - ca_vec
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rot = Bio.PDB.rotaxis(-np.pi*120.0/180.0, c_vec)
        # apply rotation to ca-n vector
        cb_at_origin = n_vec.left_multiply(rot)
        # put on top of ca atom
        cb_vec = Bio.PDB.Vector(Bio.PDB.Vector.normalized(cb_at_origin)[0:3]*CA_CB_Length[residue.get_resname()])+ca_vec 
        
        if residue.get_resname()=="GLY": # check if HA2 exits, if not, add it
            if 'HA2' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            else:
                ha2_atom = Bio.PDB.Atom.Atom('HA2',cb_vec,1,1," "," HA2",1)
                residue.add(ha2_atom)
        else:
            if 'CB' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            else:
                cb_atom = Bio.PDB.Atom.Atom('CB',cb_vec,1,1," "," CB ",1)
                residue.add(cb_atom)
    return prot 

def add_HA(prot):
    for resnum, residue in enumerate(prot.get_residues()):
        # get atom coordinates as vectors
        n_vec  = residue['N'].get_vector()
        c_vec  = residue['C'].get_vector() 
        ca_vec = residue['CA'].get_vector()
        # center at origin
        n_vec = n_vec - ca_vec
        c_vec = c_vec - ca_vec
        # find rotation matrix that rotates n 120 degrees along the ca-c vector
        rot = Bio.PDB.rotaxis(np.pi*120.0/180.0, c_vec)
        # apply rotation to ca-n vector
        ha_at_origin = n_vec.left_multiply(rot)
        # put on top of ca atom
        ha_vec = Bio.PDB.Vector(Bio.PDB.Vector.normalized(ha_at_origin)[0:3]*CA_H_Length[residue.get_resname()])+ca_vec
        
        if residue.get_resname()=="GLY": 
            if 'HA1' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            elif 'HA' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            elif 'H' in [atom.get_name() for atom in residue.get_atom()]:
                continue 
            else:
                ha_atom = Bio.PDB.Atom.Atom('H',ha_vec,1,1," "," H ",1)
                residue.add(ha_atom)
        else:
            if 'HA' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            elif 'H' in [atom.get_name() for atom in residue.get_atom()]:
                continue
            else:
                ha_atom = Bio.PDB.Atom.Atom('H',ha_vec,1,1," "," H ",1)
                residue.add(ha_atom)
    return prot

bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Program used for adding missing HN atoms to pdb files.',
        usage ='use "python %(prog)s --help or -h" for more information')
    parser.add_argument('inputpdb', help='[required] input pdb file')
    args = parser.parse_args()

    protein = Bio.PDB.PDBParser(QUIET=True).get_structure('x',args.inputpdb)
    prot = add_hn(protein)
    prot = add_CB(protein)
    prot = add_HA(protein) 

    io = Bio.PDB.PDBIO()
    io.set_structure(prot)
    io.save(bsnm(args.inputpdb) + '-wHN-wCB-wHA.pdb') 

if __name__ == '__main__':
    main()

