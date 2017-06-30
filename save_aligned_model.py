#!/usr/bin/env python
__author__  = 'Wang Zongan'

import os
import re
import sys
import string
import numpy as np
import Bio.PDB

def load_structure(pdbpath, frame_idx=0):
    return Bio.PDB.PDBParser(QUIET = True).get_structure("protein", pdbpath)[frame_idx]


def get_atoms(model, res_idx, atom_sel=['CA']):
    residues = [res for res in model.get_residues() if res.get_id()[0] == ' '] # exclude HET residues and water molecures
    residues = [residues[i] for i in res_idx]
    return [res[atom] for atom in atom_sel for res in residues]


def align_two_models(ref_model, com_model, ref_res_idx, com_res_idx, atom_sel=['CA']):
    ref_atoms = get_atoms(ref_model, ref_res_idx, atom_sel)
    com_atoms = get_atoms(com_model, com_res_idx, atom_sel)
    assert len(ref_atoms) == len(com_atoms)
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, com_atoms)
    super_imposer.apply(com_model.get_atoms())
    print "RMSD : %6.3f" % super_imposer.rms
    return com_model 


def save_model(model, outputname):
    io = Bio.PDB.PDBIO()
    io.set_structure(model) 
    io.save(outputname)
    return True


def get_fasta_length(model):
    n = 0
    for res in model.get_residues():
        if res.get_id()[0] == ' ':
            atom_nm = [atom.get_name() for atom in res.get_atom()]
            if 'CA' in atom_nm:
                n += 1
    return n 


def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
    import re

    if re.match('^([0-9]+(-[0-9]+)?)(,[0-9]+(-[0-9]+)?)*$', s) is None:
        raise ValueError('segments must be of the form 10-30,45,72-76 or similar')

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


def get_basename_without_extension(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    ref_pdb_path = sys.argv[1]
    com_pdb_path = sys.argv[2]
    ref_model = load_structure(ref_pdb_path, 0) 
    com_model = load_structure(com_pdb_path, 0)
    
    if sys.argv[3] == 'all':
        ref_res_idx = parse_segments('0-%i' % (get_fasta_length(ref_model)-1))
    else:
        ref_res_idx = parse_segments(sys.argv[3])
    if sys.argv[4] == 'all':
        com_res_idx = parse_segments('0-%i' % (get_fasta_length(com_model)-1))
    else:
        com_res_idx = parse_segments(sys.argv[4])
    print 'ref_res_idx : ', sys.argv[3]
    print 'com_res_idx : ', sys.argv[4]

    outputname = '%s.aligned.%s.pdb' % (
        get_basename_without_extension(com_pdb_path),
        get_basename_without_extension(ref_pdb_path))

    com_model = align_two_models(ref_model, com_model, ref_res_idx, com_res_idx, ['CA'])
    save_model(com_model, outputname)
    

if __name__ == '__main__':
    main()
