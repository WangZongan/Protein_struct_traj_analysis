#!/usr/bin/env python
__author__  = 'Wang Zongan'
__version__ = '2016.09.29'

import os
import re
import sys
import string
import numpy as np
import mdtraj as md
import Bio.PDB


oneletter_threeletter = dict(
    A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE',
    K='LYS', L='LEU', M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER',
    T='THR', V='VAL', W='TRP', Y='TYR')

threeletter_oneletter = dict([(v,k) for (k,v) in oneletter_threeletter.items()])

restype_order = sorted(oneletter_threeletter.values())
restype_order.append('UNH')  # 20
restype_order.append('UCO')  # 21
restype_to_index = dict((aa,i) for i,aa in enumerate(restype_order))
index_to_restype = dict((i,aa) for i,aa in enumerate(restype_order))
n_restype = 22


def bio_load(pdbpath):
    return Bio.PDB.PDBParser(QUIET=True).get_structure('protein', pdbpath)


def calc_dij_CA(model):
    nres = 0
    ref_atoms = []
    for res in model.get_residues():
        if res.get_id()[0] == ' ':  # exclude HET residues and water molecures
            atom_nm = [atom.get_name() for atom in res.get_atom()]
            if 'CA' in atom_nm:
                nres += 1
                ref_atoms.append(res['CA'])
    return np.array([ref_atoms[i] - ref_atoms[j] for i in range(nres) for j in range(i+1, nres)])


def calc_dij_CA_traj(structure):
    return np.array([calc_dij_CA(m) for m in structure])
    

def calc_dij_CA_traj_std(dij_traj):
    ''' Calculate standard deviation. '''
    return np.std(dij_traj, axis=0)


def matrix_1d_to_2d(m1d):
    '''
    m1d: 
    [(0,1), (0,2), ..., (  0,n-1),  # n-1
            (1,2), ..., (  1,n-1),  # n-2 
                   ...
                        (n-2,n-1)]  # 1

    m1 = 0 --> x = m2 - 1
    m1 = 1 --> x = (n-1) + (m2-1) - 1
    m2 = 2 --> x = (n-1)+(n-2) + (m2-2) - 1

    m1d[x] = (m1, m2) --> x = (2n - m1 - 1)*m1/2 + (m2 - m1) - 1
    '''
    n = int(np.ceil(np.sqrt(2*len(m1d))))
    m2d = np.zeros((n,n))
    for m1 in range(n):
        for m2 in range(m1+1,n):
            m2d[m1,m2] = m1d[(2*n - m1 - 1)*m1/2 + (m2 - m1) - 1]
    return m2d
         

def find_rigid_modules(dij_traj_std, cutoff, smallest_res_grp=2):
    M = np.where(dij_traj_std <= cutoff, 0, 1) 
    M = matrix_1d_to_2d(M)
    print M

    nres = M.shape[0]
    bounds_list = []
    i = 0
    while i < nres - smallest_res_grp:
        if np.nonzero(
                M[np.ix_(np.arange(i,i+smallest_res_grp),
                         np.arange(i,i+smallest_res_grp))] )[0].size == 0:
            bounds = [i, i+smallest_res_grp-1]
            j = 1
            while i+smallest_res_grp+j < nres and np.nonzero(
                                                    M[np.ix_(np.arange(i,i+smallest_res_grp+j),
                                                             np.arange(i,i+smallest_res_grp+j))]
                                                    )[0].size == 0:
                bounds = [i, i+smallest_res_grp+j-1]
                j += 1
            bounds_list.append(bounds)
            i += smallest_res_grp + j - 1
        else:
            i += 1

    return np.array(bounds_list)


def find_rigid_module_in_submatrix(mat, indices, cutoff=0.5, smallest_res_grp=2):
    subM = mat[np.ix_(indices, indices)]
    bounds_list_subM = find_rigid_modules(subM, cutoff, smallest_res_grp)+indices[0]
    return bounds_list_subM


def md_load(pdbpath): return md.load(pdbpath)

def calc_ss(traj): return md.compute_dssp(traj)

def get_consensus_ss_idx(ss):
    nmod, nres = ss.shape 
    return np.array([i for i in range(nres) if len(set(ss[:,i])) == 1]) 


def find_continuous_number(seq):
    ''' seq is in ascending order. '''
    import re
    full = np.arange(seq[0],seq[-1]+1)
    sseq = []  # string seq
    for n in full:
        if n in seq:
            sseq.append('o')
        else:
            sseq.append('_')
    csseq = ''  # continuous sseq 
    for m in re.finditer(r"o+",''.join(sseq)):
        if full[m.end()-1] > full[m.start()]:
            csseq += '%d-%d,' % (full[m.start()], full[m.end()-1])
        else:
            csseq += '%d,' % full[m.start()]
    return csseq.strip(',')


def find_nonsingular_segment(continuous_num_string):
    continuous_num_string = continuous_num_string.split(',')
    nonsingular_segments  = []
    for st in continuous_num_string:
        st = [int(_) for _ in st.split('-')]
        if len(st) > 1:
            nonsingular_segments.append(st)
    return np.array(nonsingular_segments)
    #return np.array([np.arange(a,b+1) for [a,b] in nonsingular_segments]) 


def boundary_to_range(boundaries):
    ''' boundaries = [[a0,b0],[a1,b1],...,[an-1,bn-1]] '''
    return np.concatenate(np.array([np.arange(a,b+1) for [a,b] in boundaries]))


def intersect_2boundarylist(bd1, bd2):
    inter_bd = []
    for [a1,b1] in bd1:
        for [a2,b2] in bd2:
            inter = np.intersect1d(np.arange(a1,b1+1), np.arange(a2,b2+1)) 
            if len(inter) > 1:
                inter_bd.append(find_nonsingular_segment(find_continuous_number(inter))[0])
    return inter_bd


def get_basenm_without_ext(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    pdbpath = sys.argv[1]
    cutoff  = float(sys.argv[2])

    # consensus secondary structure 
    traj = md_load(pdbpath)
    ss   = calc_ss(traj)
    con_ss_bd = find_nonsingular_segment(find_continuous_number(get_consensus_ss_idx(ss)))

    # consensus distance pairs --> local rigid fragment (LRF)
    structure = bio_load(pdbpath)
    print '%i models in trajectory.' % len(structure)
    con_loc_rig_frag_bd = find_rigid_modules(calc_dij_CA_traj_std(calc_dij_CA_traj(structure)), cutoff)

    print con_ss_bd
    print con_loc_rig_frag_bd
    with open(get_basenm_without_ext(pdbpath)+'.consensus_ss_LRF.restraint_config','w') as f:
        con_bd = intersect_2boundarylist(con_ss_bd, con_loc_rig_frag_bd) 
        if len(con_bd) == 0:
            print 'Nope, no consensus local rigid fragment (cLRF) of the same secondary structure element (sSSE) is found.'
            print >> f, ''
        else:
            print 'Consensus local rigid fragment (cLRF) of the same secondary structure element (sSSE) is found: '
            print con_bd
            output = ''
            for bd in con_bd:
                output += '--restraint-group=%s-%s ' % (bd[0],bd[1])
            print >> f, output
    

if __name__ == "__main__":
    main()


