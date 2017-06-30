# -*- encoding: utf-8 -*-
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *

import re,gc
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np

oneletter_threeletter = dict(
        A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

restype_order = sorted(oneletter_threeletter.values())
restype_order.append('UNH')  # 20
restype_order.append('UCO')  # 21

restype_to_index = dict((aa,i) for i,aa in enumerate(restype_order))
index_to_restype = dict((i,aa) for i,aa in enumerate(restype_order))
n_restype = 22
    
def seq_to_matrix(seq, seqNH, seqCO):
    nres = len(seq)
    mat = np.zeros((nres,n_restype)) 
    for i,s in enumerate(seq):
        mat[i,restype_to_index[s]] = 1.
    mat[:,20] = np.where(seqNH==0, 1, 0)
    mat[:,21] = np.where(seqCO==0, 1, 0)
    return mat
   

def load_data(data_file_path):
    Protein = collections.namedtuple(
        'Protein','nres thickness seq CAcoords CB_burial SASA'.split()) 

    colnames = 'Index chain_id chain_res_id restype secseq CB_burial NH_burial CO_burial SASA ' + \
            'NH_hbond CO_hbond SC_hbond CA_hbond CAx CAy CAz CBx CBy CBz SCMx SCMy SCMz'

    with open(data_file_path, 'r') as f:
        thickness = float(f.readline().split()[1]) 
        csv = pd.read_csv(data_file_path, delim_whitespace=True, skiprows=[0])
        ## read data ##       
        restype  = csv.restype.as_matrix()
        nres = len(restype)
        CAcoords = np.hstack((csv.CAx.as_matrix()[:,None],
                              csv.CAy.as_matrix()[:,None],
                              csv.CAz.as_matrix()[:,None]))
        NH_hbond = np.where(csv.NH_hbond.as_matrix()>0, 1, 0)
        CO_hbond = np.where(csv.CO_hbond.as_matrix()>0, 1, 0)
        CB_burial = csv.CB_burial.as_matrix()
        SASA = csv.SASA.as_matrix()
        # nres thickness seq SCMcoords include dropout CBburial  
        d = Protein(nres, 
                    thickness,
                    seq_to_matrix(restype, NH_hbond, CO_hbond),
                    CAcoords, CB_burial, SASA)
    return d


def plot_ca_projection(d,key):
    z = d.CAcoords[:,2]
    thickness = d.thickness
    idx_in_bilayer = np.where((z>=-thickness/2.)&(z<=thickness/2.))
    n = len(idx_in_bilayer[0])
    x = d.CAcoords[:,0][idx_in_bilayer]
    y = d.CAcoords[:,1][idx_in_bilayer]
    cb = d.CB_burial[idx_in_bilayer]
    sasa = d.SASA[idx_in_bilayer]

    fig, ax = subplots(figsize=(5,5), nrows=1, ncols=1)
    ax.set_title('%s : %i / %i'%(key, n, d.nres), fontsize=20)
    ax.axvline(0, color='r')
    ax.axhline(0, color='r')
    im = ax.scatter(x, y, c=cb, cmap=get_cmap('jet')) 
    lim = np.absolute([int(x.min())-5,int(x.max())+5,int(y.min())-5,int(y.max())+5]).max()
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.grid()
    fig.subplots_adjust(right=0.80)
    cbar_ax = fig.add_axes([0.80, 0.20, 0.02, 0.6]) # [left, bottom, width, height]
    cbar_ax.tick_params(labelsize=15) 
    fig.colorbar(im, cax=cbar_ax)
    tight_layout(rect=[0, 0, 0.8, 1]) # default is (0, 0, 1, 1) [left, bottom, right, top]
    savefig('TMH_testset_CAprojection_%s.png'%key,format='png')
    savefig('TMH_testset_CAprojection_%s.pdf'%key,format='pdf')


def exclude_residue_in_cavity(d, key, radius, cb_cutoff):
    z = d.CAcoords[:,2]
    thickness = d.thickness
    x = d.CAcoords[:,0]
    y = d.CAcoords[:,1]
    r = np.sqrt(x**2+y**2)
    cb = d.CB_burial
    idx_exclude = np.where(((z>=-thickness/2.)&(z<=thickness/2.)) 
                            & (r <= radius)
                            & (cb <= cb_cutoff)
                          )
    print idx_exclude
    return idx_exclude


def main():
    data_file_path = sys.argv[1]
    key = sys.argv[2]
    radius = float(sys.argv[3])
    cb_cutoff = int(sys.argv[4])

    d = load_data(data_file_path)
    plot_ca_projection(d, key)
    idx_exclude = exclude_residue_in_cavity(d, key, radius, cb_cutoff)
    with open('%s.residue_in_cavity.radius_%i.cbburialcutoff_%i.pkl'%(key, radius, cb_cutoff),'w') as f:
        cp.dump(idx_exclude[0], f, -1)

if __name__ == '__main__':
    main()
