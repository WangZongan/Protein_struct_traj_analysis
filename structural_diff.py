#!/usr/bin/env python
# coding: utf-8
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
import numpy as np
import mdtraj as md 
import Bio.PDB 

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *

__author__  = 'Wang.Zongan'
__version__ = '2016.12.14'
__description__ = '''
Given 2 pdb files with the same number of residues, which are the 2 states of the same protein, 
we want to do morphing from one state to another.

This program calculates the residue segments that have the same SS by DSSP definition 
and remain structually static under a predetermined threshold. 

We will restrain these residue segments in Upside simulation as rigid bodies. 
'''

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


def make_dij(structure, selection='CA'):
    protein = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', structure)
    print '%s Loaded.'%structure
    nres = 0
    ref_atoms = []
    for res in protein.get_residues():
        if res.get_id()[0] == ' ':  # exclude HET residues and water molecures
            atom_nm = [atom.get_name() for atom in res.get_atom()]
            if 'CA' in atom_nm:
                nres += 1
                ref_atoms.append(res[selection])
    print nres
    dist_mat = np.zeros((nres,nres))
    for i, atomi in enumerate(ref_atoms):
        for j, atomj in enumerate(ref_atoms):
            dist_mat[i][j] = atomi - atomj
    return dist_mat


def calc_diff_dist_matrix(structure1, structure2, selection='CA'):
    dist_mat1 = make_dij(structure1, selection)
    dist_mat2 = make_dij(structure2, selection)
    assert dist_mat1.shape == dist_mat2.shape
    return dist_mat1-dist_mat2


def find_rigid_modules(diff_dist_mat, cutoff, smallest_frag, inter_frag=3):
    '''
    Given a distance cutoff, calculate the modules that has internal movements no larger than the cutoff. 
    And those modules are the residue groups which are suggested to apply the restraints on. 
    '''
    M = np.where(diff_dist_mat <= cutoff, 0, 1)
    nres = diff_dist_mat.shape[0]
    bounds_list = []
    i = 0
    while i < nres - smallest_frag:
        if np.nonzero( M[np.ix_(np.arange(i, i + smallest_frag), np.arange(i, i + smallest_frag))] )[0].size == 0:
            bounds = [i, i+smallest_frag-1]
            j = 1
            while i + smallest_frag + j < nres and np.nonzero(
                        M[np.ix_(np.arange(i, i + smallest_frag + j), np.arange(i, i + smallest_frag + j))] )[0].size == 0:
                bounds = [i, i + smallest_frag + j - 1]
                j += 1
            bounds_list.append(bounds)
            i += smallest_frag + j - 1 + inter_frag 
        else:
            i += 1
    return np.array(bounds_list)


def find_rigid_module_in_diff_dist_matrix(diff_dist_mat, indices, cutoff=0.5, smallest_frag=3, inter_frag=3):
    subM = diff_dist_mat[np.ix_(indices, indices)]
    return find_rigid_modules(subM, cutoff, smallest_frag, inter_frag)+indices[0] 


def calc_ss(pdbfile):
    return md.compute_dssp(md.load(pdbfile)) 


def calc_HEC_segment_idx(pdbfile, ss='H'):
    # ss = 'H', 'E', 'C'
    return np.where((calc_ss(pdbfile)[0]==ss))[0]


def calc_same_ss(pdbfile1, pdbfile2):
    ss1 = calc_ss(pdbfile1)[0]
    ss2 = calc_ss(pdbfile2)[0]
    return np.where((ss1 == ss2))[0]


def calc_diff_ss(pdbfile1, pdbfile2):
    ss1 = calc_ss(pdbfile1)[0]
    ss2 = calc_ss(pdbfile2)[0]
    return np.where((ss1 != ss2))[0]
    

def find_continuous_number(seq):
    ''' seq is in ascending order. '''    
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


def find_nonsingular_segment(continuous_num_string):
    continuous_num_string = continuous_num_string.split(',')
    nonsingular_segments  = []
    for st in continuous_num_string:
        st = [int(_) for _ in st.split('-')]
        if len(st) > 1:
            nonsingular_segments.append(st)
    return np.array([np.arange(a,b+1) for [a,b] in nonsingular_segments])


bsnm = lambda filepath: os.path.splitext(os.path.basename(filepath))[0]


def plot_diff_dist_map(diff_dist_list, filename, distcutoff, rectanglebounds):
    from matplotlib.patches import Rectangle
    
    nres = diff_dist_list.shape[0]
    spacing = nres // 10
    def make_ticklabels(nres, spacing):
        xtl = []; ytl = []
        for i in range(nres):
            if i % spacing == 0:
                xtl.append(str(i)); ytl.append(str(i))
            else:
                xtl.append(''); ytl.append('')
        return xtl,ytl
    xtl, ytl = make_ticklabels(nres,spacing)

    fig, ax = subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(12,10))
    fig.suptitle('CA difference distance matrix', fontsize=25, x=0.435, y=0.98)
    my_cmap = get_cmap('Reds')
    my_cmap.set_under('w')
    cax = ax.pcolor(np.absolute(diff_dist_list), vmin=distcutoff, cmap=my_cmap)
    xticks(np.arange(nres), xtl)
    yticks(np.arange(nres), ytl)
    xlim(0,nres)
    ylim(0,nres)

    for bounds in rectanglebounds:
        lb, ub = bounds
        ax.add_patch(Rectangle((lb, lb), ub-lb, ub-lb, alpha=1, linewidth=3,
                    color='b', facecolor='none', fill=None))  # ower left at (x, y)
    cbar = fig.colorbar(cax)

    tick_params(axis='both', which='major', labelsize=15)
    fig.text(0.4, 0.03, 'Residue', ha='center', va='center', fontsize=25)
    fig.text(0.02, 0.5, 'Residue', ha='center', va='center', fontsize=25, rotation='vertical')
    grid()
    tight_layout(rect=[0.03, 0.05, 1, 0.95]) # default is (0, 0, 1, 1) [left, bottom, right, top]
    savefig('%s.png' % filename)



def main():
    pdbfile1      = sys.argv[1]
    pdbfile2      = sys.argv[2]
    cutoff        = float(sys.argv[3])
    smallest_frag = int(sys.argv[4])
    inter_frag    = int(sys.argv[5])
    
    print 'rigid cutoff  : %.1f' % cutoff
    print 'smallest frag : %i'   % smallest_frag
    print 'inter-frag    : %i'   % inter_frag 
    
    diff_dist_mat = calc_diff_dist_matrix(pdbfile1, pdbfile2)
    same_ss_idx   = calc_same_ss(pdbfile1, pdbfile2) # for distance file
    nres          = diff_dist_mat.shape[0]
    
    nonsingular_segments = find_nonsingular_segment(find_continuous_number(same_ss_idx))
    bounds = []
    for seg in nonsingular_segments:
        bd = find_rigid_module_in_diff_dist_matrix(diff_dist_mat, seg, cutoff=cutoff, smallest_frag=smallest_frag, inter_frag=inter_frag)
        if len(bd) != 0: 
            bounds.extend(bd)
    bounds = np.array(bounds) # for restrait groups

    filename = '%s.%s.Go.restraint_grp.rigidcutoff%.1f.smallest_frag%i.inter_frag%i.diff_dist_mat' % (
            bsnm(pdbfile1), bsnm(pdbfile2), cutoff, smallest_frag, inter_frag)
    plot_diff_dist_map(diff_dist_mat, filename, cutoff, bounds)
    
    # output 
    '''
    restraint group config : configuration for setting restraint groups in upside.
    same ss idx            : residue idx for those don't change SS  
    not restraint idx      : residue idx for those aren't restrained
    '''
    restrained = []
    with open('%s.%s.Go.restraint_grp.rigidcutoff%.1f.smallest_frag%i.inter_frag%i.config' % (
            bsnm(pdbfile1), bsnm(pdbfile2), cutoff, smallest_frag, inter_frag), 'w') as f:
        output = ''
        for bd in bounds:
            output += '--restraint-group=%s-%s ' % (bd[0],bd[1])
            restrained.extend(np.arange(bd[0],bd[1]+1))
        print >> f, output

    #with open(os.path.splitext(os.path.basename(pdbfile1))[0]+'.'+
    #          os.path.splitext(os.path.basename(pdbfile2))[0]+'.same_ss_idx.dat','w') as f:
    #    print >> f, find_continuous_number(same_ss_idx)

    # for phi & psi file
    #not_restrained = list(sorted(set(np.arange(nres))-set(restrained))) # == same_ss_not_restrained_and_diff_ss 
    #with open(os.path.splitext(os.path.basename(pdbfile1))[0]+'.'+
    #          os.path.splitext(os.path.basename(pdbfile2))[0]+'.not_restrained_idx.rigidcutoff%.2f.dat'%cutoff,'w') as f:
    #    print >> f, find_continuous_number(not_restrained)
    
    print '# restrained groups   : %i'      % len(bounds)
    print '# restrained residues : %i / %i' % (len(restrained), nres)


def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


if __name__ == "__main__":
    import time
    sta = time.time()
    main()
    print 'running time: %s' % sec_to_hr_min_sec(time.time() - sta)
    

