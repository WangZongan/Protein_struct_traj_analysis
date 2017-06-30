# coding: utf-8
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np
np.set_printoptions(precision=3,suppress=True)
import mdtraj as md 

__author__ = 'Wang.Zongan'
__version__ = '2017.06.18'

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

def seq_to_matrix(seq, seqNH, seqCO):
    nres = len(seq)
    mat = np.zeros((nres, n_restype)) 
    for i,s in enumerate(seq):
        mat[i,restype_to_index[s]] = 1.
    for i in range(nres):
        if   seqNH[i] == 0: mat[i,20] = 1
        elif seqCO[i] == 0: mat[i,21] = 1
    return mat


def calculate_ss(pdbfilename):
    prot = md.load(pdbfilename)
    ss   = md.compute_dssp(prot)
    secseq = ''.join((ele for ele in ss[0])) 
    return secseq,np.where((ss[0]=='H')) 


def secseq_to_restraint_groups_on_Hsegments(secseq):
    outputstring = ''
    ss = np.array([s for s in secseq])
    restraint = find_continuous_number(np.where(ss=='H')[0])
    for rt in restraint.split(','):
        outputstring += '--restraint-group=%s '%rt
    return outputstring


def chain_residue_idx(pdbfilename):
    prot = md.load(pdbfilename)
    top  = prot.topology
    nch  = prot.n_chains
    for i in range(nch):
        print 'chain %i' % i
        for res in top.chain(i).residues:
            print res.index, 
        print 


def find_continuous_number(seq):
    '''
    seq is in ascending order.
    '''    
    full = np.arange(seq[0],seq[-1]+1)
    
    sseq = []
    for n in full:
        if n in seq:
            sseq.append('o')
        else:
            sseq.append('_')
    
    csseq = ''
    for m in re.finditer(r"o+",''.join(sseq)):
        if full[m.end()-1] > full[m.start()]:
            csseq += '%d-%d,' % (full[m.start()], full[m.end()-1])
        else:
            csseq += '%d,' % full[m.start()]
     
    return csseq.strip(',')


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def main():
    pdb = sys.argv[1]

    secseq, secseqHidx = calculate_ss(pdb)
    
    print bsnm(pdb)
    print 'secseq:'
    print 'ss_H:\n', find_continuous_number(secseqHidx[0])
    with open('%s.Hsegments.dat' % bsnm(pdb),'w') as f:
        f.write(find_continuous_number(secseqHidx[0]))
    with open('%s.secseq.dat' % bsnm(pdb),'w') as f:
        f.write(secseq)
    with open('%s.restraint_groups.dat' % bsnm(pdb),'w') as f:
        outputstring = secseq_to_restraint_groups_on_Hsegments(secseq)
        f.write(outputstring)


if __name__ == "__main__":
    main()
