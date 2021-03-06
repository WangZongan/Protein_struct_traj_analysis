# coding: utf-8
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import numpy as np
np.set_printoptions(precision=3,suppress=True)

__author__ = 'Wang.Zongan'
__version__ = '2016.09.28'  # shortened version of protein_info.py

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


def calculate_hbond(pdbfilename):
    '''
    Geometric criteria of hbond identification given by paper: Protein Eng. 15 (2002) 359.
    D2-D1-D-H...A-A1
        Rule 1: D-A     < 3.9 A;
        Rule 2: H-A     < 2.5 A;
        Rule 3: D-H-A   > 90.0 degree;
        Rule 4: A1-A-D  > 90.0 degree;
        Rule 5: A1-A-H  > 90.0 degree.
    0,1,2 for each residue (only consider backbone-backbone hbond)
    In rare cases, nhbond = 3: bifurcated hbond
    '''
    # calculate H-bond 
    import prody
    structure = prody.parsePDB(pdbfilename)
    protein   = structure.select('protein')
    hv        = protein.getHierView()
    nres      = len([s for s in protein.select('name CA').getResnames()])

    NHbond    = np.zeros(nres,dtype=int)
    CObond    = np.zeros(nres,dtype=int)
    nhbond    = np.zeros(nres,dtype=int)
    resTyp    = ['' for i in range(nres)]

    for i, res in enumerate(hv.iterResidues()):
        resIndex         = i
        resNum           = res.getResnum()
        resTyp[resIndex] = res.getResname()
        resChainIndex    = res.getChid()
        N  = res.select('name N')
        C  = res.select('name C')
        O  = res.select('name O')
        HN = res.select('name HN')
        # N-HN
        for k,res2 in enumerate(hv.iterResidues()):
            if res2.getResindex() != resIndex:
                AO = res2.select('name O')
                if (HN != None) and (AO != None) and (prody.calcDistance(N,AO) <= 3.9):
                    if prody.calcDistance(HN,AO) <= 2.5:           # dist_H_O_check
                        if prody.calcAngle(N,HN,AO) > 90:          # angle_N_HN_O_check:
                            AC = res2.select('name C')
                            if prody.calcAngle(AC,AO,N) > 90:      # angle_C_O_N_check
                                if prody.calcAngle(AC,AO,HN) > 90: # angle_C_O_HN_check
                                    NHbond[resIndex] += 1
            if NHbond[resIndex] == 1: break  # if NHbond[resIndex] == 1 (max), break
        # C=O
        for k,res2 in enumerate(hv.iterResidues()):
            if res2.getResindex() != resIndex:
                DHN = res2.select('name HN')
                if (DHN != None) and (O != None) and (prody.calcDistance(DHN,O) <= 2.5):
                    DN = res2.select('name N')
                    if prody.calcDistance(DN,O) <= 3.9:           # dist_N_O_check
                        if prody.calcAngle(DN,DHN,O) > 90:        # angle_N_HN_O_check
                            if prody.calcAngle(C,O,DN) > 90:      # angle_C_O_N_check_
                                if prody.calcAngle(C,O,DHN) > 90: # angle_C_O_HN_check 
                                    CObond[resIndex] += 1
            if CObond[resIndex] == 2: break  # if CObond[resIndex] == 2 (max), break
    # total nhbond
    hbond = NHbond+CObond
    return NHbond, np.where((NHbond == 0)), CObond, np.where((CObond == 0)), resTyp


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


def get_basename_without_extension(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    pdb = sys.argv[1]
    NHbond, NHbond_0, CObond, CObond_0, resTyp = calculate_hbond(pdb)
    
    print 'NHbond_0 (%i residues) : \n' % len(NHbond_0[0]), find_continuous_number(NHbond_0[0])
    with open('%s.NHbond_0.dat' % get_basename_without_extension(pdb),'w') as f: 
        f.write(find_continuous_number(NHbond_0[0]))

    print 'CObond_0 (%i residues) : \n' % len(CObond_0[0]), find_continuous_number(CObond_0[0]), '\n'
    with open('%s.CObond_0.dat' % get_basename_without_extension(pdb),'w') as f:
        f.write(find_continuous_number(CObond_0[0]))
        

if __name__ == "__main__":
    main()
