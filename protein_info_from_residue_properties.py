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
__version__ = '2017.02.17'

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


#### <----- ----- ----- ----- burial cutoff ----- ----- ----- -----> ####
cb_burial_cutoff_25 = {
    'ALA': 48, 'ARG': 26, 'ASN': 34, 'ASP': 25, 'CYS': 74, 'GLN': 36, 'GLU': 24, 'GLY':  0, 'HIS': 39, 'ILE': 41,
    'LEU': 38, 'LYS': 17, 'MET': 59, 'PHE': 42, 'PRO': 27, 'SER': 42, 'THR': 43, 'TRP': 30, 'TYR': 53, 'VAL': 42}

cb_burial_cutoff_50 = {
    'ALA': 90, 'ARG': 58, 'ASN': 73, 'ASP': 60, 'CYS':103, 'GLN': 71, 'GLU': 56, 'GLY':  6, 'HIS': 74, 'ILE': 77,
    'LEU': 75, 'LYS': 44, 'MET': 92, 'PHE': 79, 'PRO': 70, 'SER': 83, 'THR': 84, 'TRP': 69, 'TYR': 88, 'VAL': 80}

cb_burial_cutoff_75 = {
    'ALA':118, 'ARG': 92, 'ASN':106, 'ASP': 93, 'CYS':124, 'GLN':100, 'GLU': 95, 'GLY': 10, 'HIS': 99, 'ILE':109,
    'LEU':106, 'LYS': 78, 'MET':116, 'PHE':107, 'PRO':107, 'SER':114, 'THR':113, 'TRP':101, 'TYR':111, 'VAL':111}

cb_burial_cutoff_100 = {
    'ALA':166, 'ARG':155, 'ASN':165, 'ASP':150, 'CYS':162, 'GLN':159, 'GLU':163, 'GLY': 22, 'HIS':152, 'ILE':163,
    'LEU':162, 'LYS':151, 'MET':159, 'PHE':156, 'PRO':168, 'SER':161, 'THR':159, 'TRP':154, 'TYR':157, 'VAL':159}

burial_cutoff = dict()
for residx in range(20):
    aa = restype_order[residx]
    burial_cutoff[aa] = np.array([cb_burial_cutoff_25[aa], cb_burial_cutoff_50[aa],
                                  cb_burial_cutoff_75[aa], cb_burial_cutoff_100[aa]])
burial_cutoff['UNH'] = np.array([11,13,15,25])
burial_cutoff['UCO'] = np.array([10,14,16,26])

def get_burial_cutoff(aa, percent):
    if   percent == 25:
        return burial_cutoff[aa][0]
    elif percent == 50:
        return burial_cutoff[aa][1]
    elif percent == 75:
        return burial_cutoff[aa][2]
    elif percent == 100:
        return burial_cutoff[aa][3]


#### <----- ----- ----- ----- hbond update ----- ----- ----- -----> ####
def assign_hbond_score_from_strength_file(csv, nres):
    hbond_scores = np.zeros((nres, 4), dtype='int') # NH, CO, CA, SC
    for i in range(len(csv)):
        if csv.H_names[i] == 'HN':
            hbond_scores[csv.donor_res_index_in_models[i]][0] += csv.relative_strengths[i]  # NH
        elif csv.H_names[i] == 'H' or csv.H_names[i] == 'HA2':
            hbond_scores[csv.donor_res_index_in_models[i]][2] += csv.relative_strengths[i]  # CA
        else:
            hbond_scores[csv.donor_res_index_in_models[i]][3] += csv.relative_strengths[i]  # SC

        if csv.A_names[i] == 'O':
            hbond_scores[csv.acceptor_res_index_in_models[i]][1] += csv.relative_strengths[i]  # CO
        else:
            hbond_scores[csv.acceptor_res_index_in_models[i]][3] += csv.relative_strengths[i]  # SC
    return hbond_scores


#### <----- ----- ----- ----- apply burial cutoff & update hbond ----- ----- ----- -----> ####
def load_d(residue_dat, hbond_dat, burialcutoff):
    with open(residue_dat, 'r') as f:
        thickness = float(f.readline().split()[1])

    csv      = pd.read_csv(residue_dat, delim_whitespace=True, skiprows=[0])
    restype  = csv.restype.as_matrix()
    nres     = len(restype)

    # hbond update
    hbond_scores = assign_hbond_score_from_strength_file(pd.read_csv(hbond_dat, delim_whitespace=True), nres)
    csv.NH_hbond = hbond_scores[:,0]
    csv.CO_hbond = hbond_scores[:,1]
    csv.CA_hbond = hbond_scores[:,2]
    csv.SC_hbond = hbond_scores[:,3]

    # burial cutoff
    is_included_sc = np.array([1 if csv.CB_burial[i] <= get_burial_cutoff(restype[i], burialcutoff) else 0 for i in range(nres)]).astype('int')
    is_included_nh = np.array([1 if csv.NH_burial[i] <= get_burial_cutoff(     'UNH', burialcutoff) else 0 for i in range(nres)]).astype('int')
    is_included_co = np.array([1 if csv.CO_burial[i] <= get_burial_cutoff(     'UCO', burialcutoff) else 0 for i in range(nres)]).astype('int')

    return restype, seq_to_matrix(restype, csv.NH_hbond, csv.CO_hbond), np.column_stack((is_included_sc, is_included_nh, is_included_co))


#### <----- ----- ----- ----- ----- ----- ----- -----> ####
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
    pdb           = sys.argv[1]
    residue_dat   = sys.argv[2]
    valid_hb_dat  = sys.argv[3]
    burial_cutoff = int(sys.argv[4])

    resTyp, seq_mat, is_included = load_d(residue_dat, valid_hb_dat, burial_cutoff)
    secseq, secseqHidx = calculate_ss(pdb)
    
    #<----- residue ----->#
    print bsnm(pdb)
    print '-'* 60
    print 'No. residues: ', len(resTyp)
    with open('%s.nres.dat' % bsnm(pdb),'w') as f:
        f.write(str(len(resTyp)))
    
    print 'resTyp:'
    ps = ''.join([threeletter_oneletter[aa] for aa in resTyp])
    for i in range(len(ps)/50+1):
        print ps[i*50:(i+1)*50]  
    
    #<----- secseq ----->#
    print 'secseq:'
    for i in range(len(ps)/50+1):
        print secseq[i*50:(i+1)*50]
    print 'ss_H:\n', find_continuous_number(secseqHidx[0])
    with open('%s.Hsegments.dat' % bsnm(pdb),'w') as f:
        f.write(find_continuous_number(secseqHidx[0]))
    with open('%s.secseq.dat' % bsnm(pdb),'w') as f:
        f.write(secseq)
    with open('%s.restraint_groups.dat' % bsnm(pdb),'w') as f:
        outputstring = secseq_to_restraint_groups_on_Hsegments(secseq)
        f.write(outputstring)

    #<----- hbond ----->#
    included_unh_idx = np.where((seq_mat[:,20]*is_included[:,1] == 1))[0]
    print 'included_unh_idx (%i): \n' % len(included_unh_idx), find_continuous_number(included_unh_idx)
    with open('%s.included_unh_idx.dat' % bsnm(pdb),'w') as f: 
        f.write(find_continuous_number(included_unh_idx))

    included_uco_idx = np.where((seq_mat[:,21]*is_included[:,2] == 1))[0]
    print 'included_uco_idx (%i): \n' % len(included_uco_idx), find_continuous_number(included_uco_idx), '\n'
    with open('%s.included_uco_idx.dat' % bsnm(pdb),'w') as f:
        f.write(find_continuous_number(included_uco_idx))

    #<----- excluded SC ----->#
    SC_excluded_idx = np.where((is_included[:,0] == 0))[0]
    print SC_excluded_idx
    if len(SC_excluded_idx) > 0:
        print 'SC excluded idx (%i): \n' % len(SC_excluded_idx), find_continuous_number(SC_excluded_idx)
        with open('%s.sc_excluded_idx.dat' % bsnm(pdb), 'w') as f:
            f.write(find_continuous_number(SC_excluded_idx))
    else:
        print 'All SCs are included.'
        with open('%s.sc_excluded_idx.dat' % bsnm(pdb), 'w') as f:
            f.write('')

        
if __name__ == "__main__":
    main()




