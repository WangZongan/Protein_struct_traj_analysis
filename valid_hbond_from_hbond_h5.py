#!/usr/bin/python
'''
Given a pdb file, up to user's choice, extract the residue properties of the protein.
'''

__author__  = 'Wang Zongan'
__version__ = '2017.02.18'

import re, gc
import sys, os, time, math, string 
import cPickle as cp
import tables  as tb 

from glob import glob
import collections
import pandas as pd
import numpy  as np
import scipy  as sp 

import Bio.PDB  # hierarchy: structure --> model --> chain --> residue --> atom 
from Bio.PDB.Vector import calc_angle 
##<----- ----- ----- ----- ----- ----- ----- ----->##

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

bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]

##<----- ----- ----- ----- read hbond h5 file ----- ----- ----- ----->##
hbond_str   = 'A1_coords A1_names A2_coords A2_names A_coords A_names D_coords D_names H_coords H_names '
hbond_str  += 'acceptor_chain_ids acceptor_res_index_in_chains acceptor_res_index_in_models acceptor_resnames '
hbond_str  += 'donor_chain_ids donor_res_index_in_chains donor_res_index_in_models donor_resnames'
hbond_tuple = collections.namedtuple('hbond', hbond_str.split())

def load_hbond_h5(fpath):
    lb = tb.open_file(fpath, mode='r')
        
    A1_coords = lb.root.A1_coords[:]
    A1_names  = lb.root.A1_names[:]
    A2_coords = lb.root.A2_coords[:]
    A2_names  = lb.root.A2_names[:]
    A_coords  = lb.root.A_coords[:]
    A_names   = lb.root.A_names[:]
    D_coords  = lb.root.D_coords[:]
    D_names   = lb.root.D_names[:]
    H_coords  = lb.root.H_coords[:]
    H_names   = lb.root.H_names[:]
                                                    
    acceptor_chain_ids           = lb.root.acceptor_chain_ids[:]
    acceptor_res_index_in_chains = lb.root.acceptor_res_index_in_chains[:]
    acceptor_res_index_in_models = lb.root.acceptor_res_index_in_models[:]
    acceptor_resnames            = lb.root.acceptor_resnames[:]
    
    donor_chain_ids           = lb.root.donor_chain_ids[:]
    donor_res_index_in_chains = lb.root.donor_res_index_in_chains[:]
    donor_res_index_in_models = lb.root.donor_res_index_in_models[:]
    donor_resnames            = lb.root.donor_resnames[:]
    
    lb.close()

    return hbond_tuple(
                A1_coords, A1_names, A2_coords, A2_names, A_coords, A_names, D_coords, D_names, H_coords, H_names,
                acceptor_chain_ids, acceptor_res_index_in_chains, acceptor_res_index_in_models, acceptor_resnames, 
                donor_chain_ids, donor_res_index_in_chains, donor_res_index_in_models, donor_resnames)

##<----- ----- ----- ----- hbond geometry ----- ----- ----- ----->##
# hbond relative strength 
calc_vec_length = lambda vec: np.sqrt(np.sum(vec*vec,axis=-1))
normalize       = lambda vec: vec/calc_vec_length(vec)

def calc_sp2_lone_pair_direction(o_atom_coord,a1_atom_coord,a2_atom_coord):    
    '''
    a2
      \     / e2 
      a1 = o 
     /     \ e1 
    '''
    norm_vec_a2_a1 = normalize(a1_atom_coord - a2_atom_coord)  # e1 
    norm_vec_a1_o  = normalize(o_atom_coord - a1_atom_coord)
    return norm_vec_a2_a1, normalize(norm_vec_a1_o-norm_vec_a2_a1)

def calc_sp3_lone_pair_direction(o_atom_coord,a1_atom_coord,oh_atom_coord):  
    '''
    rotate norm_vec_o_a1 around axis o_oh 120 degree and 240 degree to get e1 and e2
    a1    e1
      \ /
       o - e2 
       |
       oh 
    '''
    norm_vec_o_a1 = Bio.PDB.Vector.normalized(Bio.PDB.Vector(a1_atom_coord) - Bio.PDB.Vector(o_atom_coord))
    norm_vec_o_oh = Bio.PDB.Vector.normalized(Bio.PDB.Vector(oh_atom_coord) - Bio.PDB.Vector(o_atom_coord))
    rot = Bio.PDB.rotaxis(np.pi*120/180, norm_vec_o_oh)
    e1 = norm_vec_o_a1.left_multiply(rot)
    e2 = e1.left_multiply(rot) 
    return e1[:], e2[:]

d0 = 2.08
dmax = 3.9
sigmad = 0.67; sigmad2 = sigmad**2

alpha_max = 49; cos_alpha_max = np.cos(np.deg2rad(alpha_max)) 
beta_max = 37; cos_beta_max = np.cos(np.deg2rad(beta_max))

def prevent_outside_value_for_arccos(value):
    return np.select([(value>=0.99999),(value<0.99999)&(value>-0.99999),(value<=-0.99999)],
                     [0.99999,value,-0.99999])


def calc_hbond_angles(ht): 
    '''
    sidechain 
        carbonyl O:             hydroxyl O:
            ASN: OD1;           SER: OG;  
            ASP: OD1, OD2;      THR: OG1; 
            GLN: OE1;           TYR: OH.
            GLU: OE1, OE2.
    '''
    hp_strengths = np.zeros(ht.A1_coords.shape[0])
    angle_AH_lone_pair_list = np.zeros(ht.A1_coords.shape[0],dtype='f4')
    angle_AH_DH_list        = np.zeros(ht.A1_coords.shape[0],dtype='f4')
    angle_AH_AA1_list       = np.zeros(ht.A1_coords.shape[0],dtype='f4')
                        
    def calc_alpha(e_):
        return np.rad2deg(np.arccos(prevent_outside_value_for_arccos(np.vdot(n,e_))))

    for i in range(ht.A1_coords.shape[0]): 
        if ht.A2_names[i] == 'H':
            e1, e2 = calc_sp3_lone_pair_direction(
                ht.A_coords[i], ht.A1_coords[i], ht.A2_coords[i])
        else:
            e1, e2 = calc_sp2_lone_pair_direction(
                ht.A_coords[i], ht.A1_coords[i], ht.A2_coords[i])
    
        n = ht.H_coords[i] - ht.A_coords[i]  
        length_n = calc_vec_length(n)   
        n = normalize(n)
                            
        e0 = normalize(ht.H_coords[i] - ht.D_coords[i]) 
    
        norm_vec_A_A1 = normalize(ht.A1_coords[i] - ht.A_coords[i])
    
        angle_AH_lone_pair_list[i] = calc_alpha(e1) if calc_alpha(e1) <= calc_alpha(e2) else calc_alpha(e2)
        angle_AH_DH_list[i]  = np.rad2deg(np.arccos(prevent_outside_value_for_arccos(np.vdot(-n,e0))))
        angle_AH_AA1_list[i] = np.rad2deg(np.arccos(prevent_outside_value_for_arccos(np.vdot(n,norm_vec_A_A1))))
                                                                                                
    return np.array(angle_AH_lone_pair_list), np.array(angle_AH_DH_list), np.array(angle_AH_AA1_list)


def calc_hbond_position(hcoord, ocoord):
    return (1.0079*hbond+15.999*ocoord)/(1.0079+15.999)


##<----- ----- ----- ----- valid hbond ----- ----- ----- ----->##
def relative_hbond_strength(hbond_len, angle_ah_dh, angle_ah_lone):
    assert len(hbond_len) == len(angle_ah_dh) == len(angle_ah_lone)
    return (np.exp(-(E_len_func(hbond_len)-E_len_min))*
            np.exp(-(E_AH_DH_func(angle_ah_dh)-E_AH_DH_min))*
            np.exp(-(E_AH_lone_func(angle_ah_lone)-E_AH_lone_min)))

def relative_hbond_strength_valid(hbond_len, angle_ah_dh, angle_ah_lone):
    assert len(hbond_len) == len(angle_ah_dh) == len(angle_ah_lone)
    return np.where((hbond_len<=2.5)&(angle_ah_dh<=37.5)&(angle_ah_lone<=70),
                    relative_hbond_strength(hbond_len, angle_ah_dh, angle_ah_lone), 0) 

def valid_hbond_idx(hbond_len, angle_ah_dh, angle_ah_lone):
    assert len(hbond_len) == len(angle_ah_dh) == len(angle_ah_lone)
    return np.where((hbond_len<=2.5)&(angle_ah_dh<=37.5)&(angle_ah_lone<=70))

def valid_hbond(hbond_len, angle_ah_dh, angle_ah_lone):
    assert len(hbond_len) == len(angle_ah_dh) == len(angle_ah_lone)
    return np.where((hbond_len<=2.5)&(angle_ah_dh<=37.5)&(angle_ah_lone<=70), 1, 0) 


##<----- ----- ----- ----- hbond strength ----- ----- ----- ----->##
hbond_strength_tuple = collections.namedtuple('hbond_strength', 
        'type length angle_AH_lone_pair angle_AH_DH angle_AH_AA1'.split())


def write_hbond_strength_d(hbond_d):
    nh_ma  = np.where((hbond_d.H_names=='HN'),1,0)
    cah_ma = np.where((hbond_d.H_names=='H')|(hbond_d.H_names=='HA2'),1,0)
    sch_ma = np.where((hbond_d.H_names!='H')&(hbond_d.H_names!='HA2')&(hbond_d.H_names!='HN'),1,0)
    oc_ma  = np.where((hbond_d.A_names=='O'),1,0)
    osc_ma = np.where((hbond_d.A_names!='O'),1,0)
                        
    length_ = np.around(np.sqrt(np.sum((hbond_d.H_coords-hbond_d.A_coords)**2,axis=1)),3)
    alphas, betas, gamas = calc_hbond_angles(hbond_d)
                                    
    type_ = np.zeros_like(length_).astype('int')
    type_[np.nonzero( nh_ma& oc_ma)] = 1  #  NH...CO
    type_[np.nonzero( nh_ma&osc_ma)] = 2  #  NH...sidechain-O
    type_[np.nonzero(cah_ma& oc_ma)] = 3  # CAH...CO
    type_[np.nonzero(cah_ma&osc_ma)] = 4  # CAH...sidehain-O
    type_[np.nonzero(sch_ma& oc_ma)] = 5  # sidechain-H...CO
    type_[np.nonzero(sch_ma&osc_ma)] = 6  # sidechain-H...sidechain-O
    
    return hbond_strength_tuple(
            type_, 
            length_, 
            alphas, betas, gamas,  # (angle_AH_lone_pair_list), (angle_AH_DH_list), (angle_AH_AA1_list)
        )


def write_valid_hbond(fpath, hbond_d):
    with open(fpath, 'w') as f:
        print >> f, 'H_names donor_res_index_in_models A_names acceptor_res_index_in_models length angle_AH_lone_pair angle_AH_DH angle_AH_AA1 relative_strengths'

        hbond_strength_d =  write_hbond_strength_d(hbond_d)
                
        valids = valid_hbond(
                    hbond_strength_d.length, 
                    hbond_strength_d.angle_AH_DH, 
                    hbond_strength_d.angle_AH_lone_pair).astype('int')

        for i in range(len(valids)):
            if valids[i] > 0:
                print >> f, '%4s %4i %4s %4i %6.3f %6.3f %6.3f %6.3f %i' % (
                    hbond_d.H_names[i], 
                    hbond_d.donor_res_index_in_models[i], 
                    hbond_d.A_names[i], 
                    hbond_d.acceptor_res_index_in_models[i],
                    hbond_strength_d.length[i], 
                    hbond_strength_d.angle_AH_lone_pair[i], 
                    hbond_strength_d.angle_AH_DH[i], 
                    hbond_strength_d.angle_AH_AA1[i], 
                    valids[i])
    return True

##<----- ----- ----- ----- ----- ----- ----- ----->##
def main():
    hbond_h5  = sys.argv[1]
    valid_dat = sys.argv[2]

    hbond_d   = load_hbond_h5(hbond_h5)
    write_valid_hbond(valid_dat, hbond_d)



##<----- ----- ----- ----- ----- ----- ----- ----->##
if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    end_time   = time.time()
    print
    print 'Total running time:'
    print "--- %s seconds ---" % round((end_time - start_time),3)
