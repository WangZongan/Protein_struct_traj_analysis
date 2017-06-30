#!/usr/bin/env python
'''
Given a pdb file, calculate the ANM from it.

'''
__author__  = 'Wang Zongan'
__version__ = '2016-11-16'

import os
import sys
import string
import numpy as np
import prody as prd

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def calc_ANM_CA_prody(CAcoords):
    anm = prd.ANM()
    anm.buildHessian(CAcoords, cutoff=15.0) # default = 15
    anm.calcModes(n_modes=None)
    return anm.getEigvals().round(6), anm.getEigvecs().round(6)
    

def calc_cumul_corr_cos(dX, anmEigvecs):
    '''
    dX         : ( N, 3)
    anmEigvecs : (3N, 3N-6) 
    '''
    return np.cumsum(np.inner(np.hstack(dX), anmEigvecs.T)**2)/np.sum(dX**2)


def rmsd_transform(target, model):
    assert target.shape == model.shape == (model.shape[0],3)
    base_shift_target = target.mean(axis=0)
    base_shift_model  = model .mean(axis=0)

    target = target - target.mean(axis=0)
    model  = model  - model .mean(axis=0)

    R = np.dot(target.T, model)
    U,S,Vt = np.linalg.svd(R)
    if np.linalg.det(np.dot(U,Vt))<0.:
        Vt[:,-1] *= -1.  # fix improper rotation
    rot = np.dot(U,Vt)
    shift = base_shift_target - np.dot(rot, base_shift_model)
    return rot, shift

def structure_rmsd(a,b):
    rot,trans = rmsd_transform(a,b)
    diff = a - (trans+np.dot(b,rot.T))
    return np.sqrt((diff**2).sum(axis=-1).mean(axis=-1))

def traj_rmsd(traj, native):
    return np.array([structure_rmsd(x,native) for x in traj])

def structure_align(a,b):
    rot,trans = rmsd_transform(a,b)
    return trans+np.dot(b,rot.T)

def traj_align(traj, native):
    return np.array([structure_align(x, native) for x in traj])


def plot_line(Ys, name_lists, fig_name, plot_format='png', plot_color_map='jet'): 
    nY = len(Ys)

    fig = figure(figsize=(10,10))
    title('Cumul Corr Cos',fontsize=20)
    my_cmap = get_cmap(plot_color_map)(np.linspace(0,1,nY))

    for i,c in zip(range(nY),my_cmap):
        plot(Ys[i], color=c, label='image %i'%name_lists[i])

    legend(loc='best',fontsize='medium')  
    xlabel('# modes',fontsize=20) 
    ylabel('Cumulative Squred Cosine',fontsize=20)
    xlim(0,)
    ylim(0,)
    tick_params(axis='both', which='major', labelsize=15)
    grid()
    tight_layout() # default is (0, 0, 1, 1) [left, bottom, right, top]
    savefig('%s.cumul_corr_cos.%s'%(fig_name, plot_format), format=plot_format) 


def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        usage ='use "python %(prog)s --help" for more information',
        description = 'Calculate the cumulative correlation cosine. ')
    parser.add_argument('pdb', help='[required] input PDB file.')
    args = parser.parse_args()

    p_ca = prd.parsePDB(args.pdb).select('protein and name CA')
    coords = p_ca.getCoordsets(None)
    eigenvals, eigenvecs = calc_ANM_CA_prody(coords[0])
    indices = np.linspace(0, coords.shape[0]-1, 10).astype(int)[1:]

    CCCs = []
    for coord in coords[indices]:
        CCCs.append(calc_cumul_corr_cos(structure_align(coord,coords[0])-coords[0], eigenvecs))

    plot_line(CCCs, indices, os.path.splitext(os.path.basename(args.pdb))[0]) 
    

if __name__ == '__main__':
    import time
    sta = time.time()
    main()
    print 'running time: %s' % sec_to_hr_min_sec(time.time() - sta)

