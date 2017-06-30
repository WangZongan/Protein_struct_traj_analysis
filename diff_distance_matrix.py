#!/usr/bin/env python
import os
import sys
import string
import numpy as np
import Bio.PDB 
import cPickle as cp

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *

def make_dij(structure, selection='CA'):
    protein = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', structure)
    print '%s Loaded.'%structure
    nres = 0
    ref_atoms = []
    for res in protein.get_residues():
        if res.get_id()[0] == ' ':  # exclude HET residues and water molecures
            nres += 1
            ref_atoms.append(res[selection])
    dist_list = np.zeros((nres,nres))
    for i, atomi in enumerate(ref_atoms):
        for j, atomj in enumerate(ref_atoms):
            dist_list[i][j] = atomi - atomj
    return dist_list 

def calc_diff_dist_matrix(structure1, structure2, selection='CA'):
    dist_list1 = make_dij(structure1, selection)
    dist_list2 = make_dij(structure2, selection)
    assert dist_list1.shape == dist_list2.shape
    return dist_list1-dist_list2

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
    savefig('%s.diff_dist_map.distcutoff%i.png'%(filename,distcutoff))


def find_rigid_modules(diff_dist_list, distcutoff, smallest_res_grp=10):
    '''
    Given a distance cutoff, calculate the modules that has internal movements no larger than the cutoff. 
    And those modules are the residue groups which are suggested to apply the restraints on. 
    '''
    M = np.where(diff_dist_list <= distcutoff, 0, 1)
    nres = diff_dist_list.shape[0]
    bounds_list = []
    i = 0
    while i < nres - smallest_res_grp:
        if np.nonzero(M[i:i+smallest_res_grp,i:i+smallest_res_grp])[0].size == 0:
            bounds = [i, i+smallest_res_grp]
            j = 1
            while i+smallest_res_grp+j < nres and np.nonzero(M[i:i+smallest_res_grp+j,i:i+smallest_res_grp+j])[0].size == 0:
                bounds = [i, i+smallest_res_grp+j]
                j += 1
            bounds_list.append(bounds)
            i += smallest_res_grp+j
        else:
            i += 1
    return np.array(bounds_list)


def main():
    pdbfile1 = sys.argv[1]
    pdbfile2 = sys.argv[2]
    cutoff = int(sys.argv[3])

    diff_dist_list = calc_diff_dist_matrix(pdbfile1, pdbfile2)
    filename = os.path.basename(pdbfile1)+'.'+os.path.basename(pdbfile2)
    bounds = find_rigid_modules(diff_dist_list, cutoff)
    plot_diff_dist_map(diff_dist_list, filename, cutoff, bounds)

    with open(filename+'.rigid_module.cutoff%i.pkl'%cutoff,'w') as f:
        cp.dump(bounds,f,-1)

    print bounds


if __name__=="__main__":
    main()

