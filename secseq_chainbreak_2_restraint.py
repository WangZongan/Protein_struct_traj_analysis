# coding: utf-8
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np

np.set_printoptions(precision=3,suppress=True)


# In[11]:
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

def seq_to_matrix(seq,seqNH,seqCO):
    nres = len(seq)
    mat = np.zeros((nres,n_restype)) 
    for i,s in enumerate(seq):
        mat[i,restype_to_index[s]] = 1.
    for i in range(nres):
        if   seqNH[i] == 0: mat[i,20] = 1
        elif seqCO[i] == 0: mat[i,21] = 1
    return mat


# In[9]:
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


# In[19]:

def main():
    secseq_file = sys.argv[1] # sec. sequence  
    chain_breaks_file = sys.argv[2] # n-1, n, n+1 (n is the first residue of a chain) 

    with open(secseq_file,'r') as f:
        secseq = f.readline()
    nres = len(secseq)

    with open(chain_breaks_file,'r') as f:
        chain_break = f.readline().split(',')
        for i in range(len(chain_break)):
            chain_break[i] = int(chain_break[i])
        chain_break = np.array(list(set(chain_break)))
        chain_break = chain_break[1::3]  # first residue in each chain

    chain_secseq = [secseq[0:chain_break[0]]]
    for i in range(len(chain_break)-1):
        chain_secseq.append(secseq[chain_break[i]:chain_break[i+1]])
    chain_secseq.append(secseq[chain_break[-1]:])

    outputstring = ''
    for i, cs in enumerate(chain_secseq):
        cs = np.array([l for l in cs])
        print np.where(cs=='H')[0]
        if i == 0:
            restraint = find_continuous_number(np.where(cs=='H')[0])
        else:
            restraint = find_continuous_number(np.where(cs=='H')[0]+chain_break[i-1])
        print restraint
        for rt in restraint.split(','):
            outputstring += '--restraint-group=%s '%rt
    
    print outputstring
    with open(secseq_file[0:4]+'.Hsegment_restraints_config.dat','w') as f:
        f.write(outputstring)
    


if __name__ == "__main__":
    main()

