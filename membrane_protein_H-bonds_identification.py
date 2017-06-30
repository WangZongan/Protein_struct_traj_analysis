#!/usr/bin/python
'''
This program is for identifing H-bonds using mdtraj from pdb file with side chains completed by psfgen on VMD.

Please be noted: the side chains should be completed by psfgen !
'''

import matplotlib
from matplotlib.pyplot import *
import seaborn as sns
sns.set_style(style='white')
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

import re,gc
import sys,os,time,math,string
import cPickle as cp
import tables as tb

from glob import glob
import collections

import pandas as pd
import numpy as np
import scipy as sp

np.set_printoptions(precision=3,suppress=True)

import mdtraj as md

sys.path.insert(0, '/Users/wangzongan/Documents/experiment/upside_latest/jobs/scripts/py_files')
import mdtraj_hbond_WZ as mdwz

__author__ = 'Wang Zongan'
__version__ = '2017.06.28'

'''
ALA
N H CA HA CB HB1 HB2 HB3 C O
ARG
N H H2 H3 CA HA CB HB3 HB2 CG HG3 HG2 CD HD3 HD2 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O
ASN
N H CA HA CB HB3 HB2 CG OD1 ND2 HD21 HD22 C O
ASP
N H CA HA CB HB3 HB2 CG OD1 OD2 C O
CYS
N H CA HA CB HB3 HB2 SG C O
GLN
N H CA HA CB HB3 HB2 CG HG3 HG2 CD OE1 NE2 HE21 HE22 C O
GLU
N H CA HA CB HB3 HB2 CG HG3 HG2 CD OE1 OE2 C O
GLY
N H CA HA3 HA2 C O
HIS
N H CA HA CB HB3 HB2 ND1 HD1 CG CE1 HE1 NE2 CD2 HD2 C O
ILE
N H CA HA CB HB CG2 HG21 HG22 HG23 CG1 HG13 HG12 CD1 HD11 HD12 HD13 C O
LEU
N H CA HA CB HB3 HB2 CG HG CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23 C O
LYS
N H CA HA CB HB3 HB2 CG HG3 HG2 CD HD3 HD2 CE HE3 HE2 NZ HZ1 HZ2 HZ3 C O
MET
N H CA HA CB HB3 HB2 CG HG3 HG2 SD CE HE1 HE2 HE3 C O
PHE
N H CA HA CB HB3 HB2 CG CD1 HD1 CE1 HE1 CZ HZ CD2 HD2 CE2 HE2 C O
PRO
N CD HD3 HD2 CA HA CB HB3 HB2 CG HG3 HG2 C O
SER
N H CA HA CB HB3 HB2 OG HG C O
THR
N H CA HA CB HB OG1 HG1 CG2 HG21 HG22 HG23 C O
TRP
N H CA HA CB HB3 HB2 CG CD1 HD1 NE1 HE1 CE2 CD2 CE3 HE3 CZ3 HZ3 CZ2 HZ2 CH2 HH2 C O
TYR
N H CA HA CB HB3 HB2 CG CD1 HD1 CE1 HE1 CZ OH HH CD2 HD2 CE2 HE2 C O
VAL
N H CA HA CB HB CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23 C O
'''

oneletter_threeletter = dict(
        A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

restype_order = sorted(oneletter_threeletter.values())
restype_order.append('UNH')  # 20
restype_order.append('UCO')  # 21
restype_order.append('UCA')  # 22
restype_to_index = dict((aa,i) for i,aa in enumerate(restype_order))
index_to_restype = dict((i,aa) for i,aa in enumerate(restype_order))
n_restype = 23

def seq_to_matrix(seq, seqNH, seqCO, seqCA):
    nres = len(seq)
    mat = np.zeros((nres, n_restype))
    for i,s in enumerate(seq):
        mat[i,restype_to_index[s]] = 1.
    for i in range(nres):
        if   seqNH[i] == 0: mat[i,20] = 1
        elif seqCO[i] == 0: mat[i,21] = 1
        elif seqCA[i] == 0: mat[i,22] = 1
    return mat


# In[128]:
hbond_donor_related_atoms_dict = {
    'ALA':[['H','N'],['HA','CA']],
    'ARG':[['H','N'],['HA','CA'],['HE','NE'],['HH11','NH1'],['HH12','NH1'],['HH21','NH2'],['HH22','NH2']],
    'ASN':[['H','N'],['HA','CA'],['HD21','ND2'],['HD22','ND2']],
    'ASP':[['H','N'],['HA','CA']],
    'CYS':[['H','N'],['HA','CA']],
    'GLN':[['H','N'],['HA','CA'],['HE21','NE2'],['HE22','NE2']],
    'GLU':[['H','N'],['HA','CA']],
    'GLY':[['H','N'],['HA2','CA'],['HA3','CA']],
    'HIS':[['H','N'],['HA','CA'],['HD1','ND1']],
    'ILE':[['H','N'],['HA','CA']],
    'LEU':[['H','N'],['HA','CA']],
    'LYS':[['H','N'],['HA','CA'],['HZ1','NZ'],['HZ2','NZ'],['HZ3','NZ']],
    'MET':[['H','N'],['HA','CA']],
    'PHE':[['H','N'],['HA','CA']],
    'PRO':[          ['HA','CA']],
    'SER':[['H','N'],['HA','CA'],['HG','OG']],
    'THR':[['H','N'],['HA','CA'],['HG1','OG1']],
    'TRP':[['H','N'],['HA','CA'],['HE1','NE1']],
    'TYR':[['H','N'],['HA','CA'],['HH','OH']],
    'VAL':[['H','N'],['HA','CA']]
}

hbond_acceptor_related_atoms_dict = {
    'ALA':[['O','C','CA'],                                       ['N','CA','H']],
    'ARG':[['O','C','CA'],                                       ['N','CA','H'], ['NE','CD','CZ'], ['NH1','HH11','CZ'], ['NH2','HH21','CZ']],
    'ASN':[['O','C','CA'], ['OD1','CG','CB'],                    ['N','CA','H'], ['ND2','CG','HD21']],
    'ASP':[['O','C','CA'], ['OD1','CG','CB'], ['OD2','CG','CB'], ['N','CA','H']],
    'CYS':[['O','C','CA'],                                       ['N','CA','H']],
    'GLN':[['O','C','CA'], ['OE1','CD','CG'],                    ['N','CA','H'], ['NE2','CD','HE21']],
    'GLU':[['O','C','CA'], ['OE1','CD','CG'], ['OE2','CD','CG'], ['N','CA','H']],
    'GLY':[['O','C','CA'],                                       ['N','CA','H']],
    'HIS':[['O','C','CA'],                                       ['N','CA','H'], ['ND1','CG','CE1'], ['NE2','CD2','CE1']],
    'ILE':[['O','C','CA'],                                       ['N','CA','H']],
    'LEU':[['O','C','CA'],                                       ['N','CA','H']],
    'LYS':[['O','C','CA'],                                       ['N','CA','H'], ['NZ','CE','HZ1']],
    'MET':[['O','C','CA'],                                       ['N','CA','H']],
    'PHE':[['O','C','CA'],                                       ['N','CA','H']],
    'PRO':[['O','C','CA'],                                       ['N','CA','CD']],
    'SER':[['O','C','CA'], ['OG','CB','HG'],                     ['N','CA','H']],
    'THR':[['O','C','CA'], ['OG1','CB','HG1'],                   ['N','CA','H']],
    'TRP':[['O','C','CA'],                                       ['N','CA','H'], ['NE1','CD1','CE2']],
    'TYR':[['O','C','CA'], ['OH','CZ','HH'],                     ['N','CA','H']],
    'VAL':[['O','C','CA'],                                       ['N','CA','H']]
}


def find_atom_info(mdtraj_top, atom_index_in_model):

    nch  = mdtraj_top.n_chains
    nres = mdtraj_top.n_residues
    nres_eachchain = [ch.n_residues for ch in mdtraj_top.chains]

    at = mdtraj_top.atom(atom_index_in_model)

    i = 0
    res_index_in_chain = at.residue.index
    while res_index_in_chain > nres_eachchain[i] and i < nch:
        res_index_in_chain -= nres_eachchain[i]
        i += 1

    return at.residue.chain.index, res_index_in_chain, at.residue.index, at.residue.name


def find_acceptor_related_atoms(mdtraj_top, A_id):

    A_at = mdtraj_top.atom(A_id)
    #print A_at.name, A_at.residue.name,
    #for iat in range(A_at.residue.n_atoms):
    #    print A_at.residue.atom(iat).name,
    #print

    def name_to_id(A1_at_nm, A2_at_nm):
        for iat in range(A_at.residue.n_atoms):
            if A_at.residue.atom(iat).name == A1_at_nm:
                A1_id = A_at.residue.atom(iat).index
            if A_at.residue.atom(iat).name == A2_at_nm:
                A2_id = A_at.residue.atom(iat).index
        return A1_id, A2_id

    if A_at.name == 'O' or A_at.name == 'OXT':
        A1_id, A2_id = name_to_id('C', 'CA')
    elif A_at.name == 'OD1' or A_at.name == 'OD2':
        A1_id, A2_id = name_to_id('CG', 'CB')
    elif A_at.name == 'OE1' or A_at.name == 'OE2':
        A1_id, A2_id = name_to_id('CD', 'CG')
    elif A_at.name == 'OG':
        A1_id, A2_id = name_to_id('CB', 'HG')
    elif A_at.name == 'OG1':
        A1_id, A2_id = name_to_id('CB', 'HG1')
    elif A_at.name == 'OH':
        A1_id, A2_id = name_to_id('CZ', 'HH')

    elif A_at.name == 'NH1':
        A1_id, A2_id = name_to_id('HH11', 'CZ')
    elif A_at.name == 'NH2':
        A1_id, A2_id = name_to_id('HH21', 'CZ')

    elif A_at.name == 'ND1':
        A1_id, A2_id = name_to_id('CG', 'CE1')
    elif A_at.name == 'ND2':
        A1_id, A2_id = name_to_id('CG', 'HD21')

    elif A_at.name == 'NE':
        A1_id, A2_id = name_to_id('CD', 'CZ')
    elif A_at.name == 'NE1':
        A1_id, A2_id = name_to_id('CD1', 'CE2')
    elif A_at.name == 'NE2' and A_at.residue.name == 'GLN':
        A1_id, A2_id = name_to_id('CD', 'HE21')
    elif A_at.name == 'NE2' and A_at.residue.name == 'HIS':
        A1_id, A2_id = name_to_id('CE1', 'CD2')

    elif A_at.name == 'NZ':
        A1_id, A2_id = name_to_id('CE', 'HZ1')

    return A1_id, A2_id


donor_info_tuple                = collections.namedtuple('donor'                    , 'chain_id res_index_in_chain res_index_in_model resname donor_atom_name'.split())
donor_related_atom_names_tuple  = collections.namedtuple('donor_related_atom_names' , 'D H'.split())
donor_related_atom_coords_tuple = collections.namedtuple('donor_related_atom_coords', 'D H'.split())

# A2 is connected to either A1 or A
acceptor_info_tuple                = collections.namedtuple('acceptor'                    ,'chain_id res_index_in_chain res_index_in_model resname acceptor_atom_name'.split())
acceptor_related_atom_names_tuple  = collections.namedtuple('acceptor_related_atom_names' , 'A A1 A2'.split())
acceptor_related_atom_coords_tuple = collections.namedtuple('acceptor_related_atom_coords', 'A A1 A2'.split())


donor_info_string    = 'donor_chain_id    donor_res_index_in_chain    donor_res_index_in_model    donor_resname '
acceptor_info_string = 'acceptor_chain_id acceptor_res_index_in_chain acceptor_res_index_in_model acceptor_resname '

hbond_tuple_string   = donor_info_string    + 'D_name H_name D_coord H_coord '
hbond_tuple_string  += acceptor_info_string + 'A_name A1_name A2_name A_coord A1_coord A2_coord'
hbond_tuple          = collections.namedtuple('HBond', hbond_tuple_string.split())


def find_hbond_in_model(mdtraj_protein,
                        hbond_identifier='baker_hubbard',
                        distance_cutoff=0.30, angle_cutoff=2.094395):
    '''
    hbond_identifier:
    baker_hubbard, wernet_nilsson, vmd_hbond

    hbonds_DHA_ids: np.array, (n_pairs, 3) (D_id, H_id, A_ id)

    At this moment, ignore CAH bonds.
    '''

    if hbond_identifier == 'baker_hubbard':
        hbonds_DHA_ids = mdwz.baker_hubbard(mdtraj_protein,
                            distance_cutoff=0.30, angle_cutoff=2.094395, periodic=False)
    elif hbond_identifier == 'vmd_hbond':
        hbonds_DHA_ids = mdwz.vmd_hbond(mdtraj_protein,
                            distance_cutoff=0.40, angle_cutoff=2.094395, periodic=False)
    elif hbond_identifier == 'wernet_nilsson':
        hbonds_DHA_ids = md.wernet_nilsson(mdtraj_protein, periodic=False)[0]

    top    = mdtraj_protein.topology
    nres   = top.n_residues
    coords = mdtraj_protein.xyz[0]

    # simple counts just to see to be or not to be
    simple_NH_hbond_count = np.zeros(nres,dtype=int)
    simple_CO_hbond_count = np.zeros(nres,dtype=int)
    simple_SC_hbond_count = np.zeros(nres,dtype=int)
    simple_CA_hbond_count = np.zeros(nres,dtype=int)

    hbonds_pairs_in_model = []

    for hbond_pair_ids in hbonds_DHA_ids:

        D_id = hbond_pair_ids[0]
        H_id = hbond_pair_ids[1]
        A_id = hbond_pair_ids[2]

        donor_ch_id, donor_res_id_ch, donor_res_id_model, donor_resnm             = find_atom_info(top, H_id)
        acceptor_ch_id, acceptor_res_id_ch, acceptor_res_id_model, acceptor_resnm = find_atom_info(top, A_id)

        if top.atom(A_id).name == 'N':
            continue
        else:
            A1_id, A2_id = find_acceptor_related_atoms(top, A_id)

        # donor
        if top.atom(D_id).name == 'N':
            simple_NH_hbond_count[donor_res_id_model] += 1
        elif top.atom(D_id).name == 'CA':
            simple_CA_hbond_count[donor_res_id_model] += 1
        else:
            simple_SC_hbond_count[donor_res_id_model] += 1

        # acceptor
        if top.atom(A_id).name == 'O':
            simple_CO_hbond_count[acceptor_res_id_model] += 1
        else:
            simple_SC_hbond_count[acceptor_res_id_model] += 1

        hbonds_pairs_in_model.append(
            hbond_tuple(
                donor_ch_id, donor_res_id_ch, donor_res_id_model, donor_resnm,  # donor_info
                top.atom(D_id).name, top.atom(H_id).name,                       # donor_related_atom_names
                coords[D_id],             # D_coord
                coords[H_id],             # H_coord
                acceptor_ch_id, acceptor_res_id_ch, acceptor_res_id_model, acceptor_resnm,  # acceptor_info
                top.atom(A_id).name, top.atom(A1_id).name, top.atom(A2_id).name,            # acceptor_related_atom_names
                coords[A_id],             # A_coord
                coords[A1_id],            # A1_coord
                coords[A2_id]             # A2_coord
            )
        )

    return hbonds_pairs_in_model, np.array([simple_NH_hbond_count, simple_CO_hbond_count,
                                            simple_SC_hbond_count, simple_CA_hbond_count])


def write_N_O_coords(mdtraj_protein, nm):
    def name_to_id(res, at_name):
        for iat in range(res.n_atoms):
            if res.atom(iat).name == at_name:
                at_id = res.atom(iat).index
        return at_id

    top    = mdtraj_protein.topology
    coords = mdtraj_protein.xyz[0]
    nres   = top.n_residues

    Ncoords = np.zeros((nres,3))
    Ocoords = np.zeros((nres,3))
    for i in range(nres):
        res = top.residue(i)
        N_id = name_to_id(res, 'N'); Ncoords[i] = coords[N_id]
        O_id = name_to_id(res, 'O'); Ocoords[i] = coords[O_id]

    with open('%s.NOcoords.pkl' % (nm), 'w') as f:
         cp.dump((Ncoords, Ocoords), f, -1)

    return True


def main():
    pdb_fpath = sys.argv[1]
    output_nm = sys.argv[2]
    hbond_identifier = sys.argv[3]
    distance_cutoff  = float(sys.argv[4]) # nm
    angle_cutoff     = float(sys.argv[5]) # degree

    test_protein = md.load_pdb(pdb_fpath)
    hbonds_pairs_in_model, hbond_counts = find_hbond_in_model(test_protein,
                                                hbond_identifier=hbond_identifier,
                                                distance_cutoff=distance_cutoff,
                                                angle_cutoff=angle_cutoff*np.pi/180)

    with tb.open_file('%s.hbond_pairs.%s.dist_cf%.2f.angle_cf%.1f.h5' %
                        (output_nm, hbond_identifier, distance_cutoff, angle_cutoff), mode='w') as lib:

        lb.create_array(lb.root, 'donor_chain_ids',              np.array([p.donor_chain_id              for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'donor_res_index_in_chains',    np.array([p.donor_res_index_in_chain    for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'donor_res_index_in_models',    np.array([p.donor_res_index_in_model    for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'donor_resnames',               np.array([p.donor_resname               for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'D_names',                      np.array([p.D_name                      for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'H_names',                      np.array([p.H_name                      for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'D_coords',                     np.array([p.D_coord                     for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'H_coords',                     np.array([p.H_coord                     for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'acceptor_chain_ids',           np.array([p.acceptor_chain_id           for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'acceptor_res_index_in_chains', np.array([p.acceptor_res_index_in_chain for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'acceptor_res_index_in_models', np.array([p.acceptor_res_index_in_model for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'acceptor_resnames',            np.array([p.acceptor_resname            for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A_names',                      np.array([p.A_name                      for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A1_names',                     np.array([p.A1_name                     for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A2_names',                     np.array([p.A2_name                     for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A_coords',                     np.array([p.A_coord                     for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A1_coords',                    np.array([p.A1_coord                    for p in hbonds_pairs_in_model]))
        lb.create_array(lb.root, 'A2_coords',                    np.array([p.A2_coord                    for p in hbonds_pairs_in_model]))

    with open('%s.hbond_counts.%s.dist_cf%.2f.angle_cf%.1f.pkl' %
                (output_nm, hbond_identifier, distance_cutoff, angle_cutoff), 'w') as f:
        cp.dump(hbond_counts, f, -1)

    write_N_O_coords(test_protein, output_nm)


if __name__ == '__main__':
    main()
