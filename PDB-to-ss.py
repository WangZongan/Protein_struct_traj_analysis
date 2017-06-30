#!/usr/bin/env python
# coding: utf-8
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
import numpy as np
import mdtraj as md 
import Bio.PDB 

__author__ = 'Wang.Zongan'
__version__ = '2016.09.26'

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

def calculate_ss(pdbfilename, simplified=True):
    '''
    The DSSP assignment codes are:
        H : Alpha helix
        B : Residue in isolated beta-bridge
        E : Extended strand, participates in beta ladder
        G : 3-helix (3/10 helix)
        I : 5 helix (pi helix)
        T : hydrogen bonded turn
        S : bend
          : Loops and irregular elements

    There are two ways to simplify 8-letter DSSP codes. 
    By default, the simplified DSSP codes in mdtraj are:
        H : Helix. Either of the H, G, or I codes.
        E : Strand. Either of the E, or B codes.
        C : Coil. Either of the T, S or ' ' codes.

    Simplify DSSP codes in this way:
        H : H
        E : E
        C : all the others
    '''
    import mdtraj as md
    prot = md.load(pdbfilename)
    ss   = md.compute_dssp(prot, simplified=False)[0]
    
    if simplified == True:
        ss[np.where((ss!='H')&(ss!='E'))] = 'C'

    return ss


def get_basename_without_ext(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    pdbfile = sys.argv[1]
    with open(get_basename_without_ext(pdbfile)+'.ss3','w') as f:
        print >> f, '> %s' % pdbfile
        print >> f, ''.join(calculate_ss(pdbfile))


if __name__ == "__main__":
    main()
