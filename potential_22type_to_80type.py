# coding: utf-8
import re,gc
import sys,os,time,math,string 
import cPickle as cp
import collections
import numpy as np
from matplotlib.pyplot import *
import matplotlib.cm as cm
import seaborn as sns
sns.set_style(style='white')
import tables as tb

__version__ = '2017.02.17'
__author__  = 'Wang Zongan'

#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----# 
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

def seq_to_matrix(seq,seqNH,seqCO):
    nres = len(seq)
    mat = np.zeros((nres,n_restype))
    for i,s in enumerate(seq):
        mat[i,restype_to_index[s]] = 1.
    mat[:,20] = np.where(seqNH==0, 1, 0) # UNH
    mat[:,21] = np.where(seqCO==0, 1, 0) # UCO
    return mat


#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----#
def read_potential_h5_file_22type(h5file):
    lib      = tb.open_file(h5file, mode='r')
    names    = lib.root.names[:]
    z_energy = lib.root.z_energy[:]
    z_min    = lib.root.z_energy._v_attrs.z_min
    z_max    = lib.root.z_energy._v_attrs.z_max
    try:
        thickness = lib.root.z_energy._v_attrs.thickness
    except AttributeError:
        thickness = 30.0
    lib.close()
    return names, z_energy, z_min, z_max, thickness


import scipy.interpolate
def extrapolated_spline_1D(x0,y0,k=1):
    spline = scipy.interpolate.InterpolatedUnivariateSpline(x0,y0,k=k)
    def f(x, spline=spline):
        return np.select(
                    [(x<x0[0]),              (x>x0[-1]),              np.ones_like(x,dtype='bool')],
                    [np.zeros_like(x)+y0[0], np.zeros_like(x)+y0[-1], spline(x)])
    return f


#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----#
def write_potential_h5_file_80type(names, z_energy, z_min, z_max, thickness, h5filename):
    forcezeros = np.array([extrapolated_spline_1D(np.linspace(z_min, z_max, len(z_energy[0])), z_energy[i], 3)(thickness/2.+20) for i in range(n_restype)])
    for i in range(n_restype):
        z_energy[i] = z_energy[i] - forcezeros[i]

    print names, thickness
    new_names = []
    for nm in names[:20]:
        new_names.extend([nm, nm+'_UNH', nm+'_UCO', nm+'_UNH_UCO'])

    new_z_energy         = np.zeros((80, len(z_energy[0])))
    new_z_energy[0:80:4] = z_energy[0:20]
    if 'UNH' in names and 'UCO' in names: 
        new_z_energy[1:80:4] = z_energy[0:20] + z_energy[20] # _UNH
        new_z_energy[2:80:4] = z_energy[0:20] + z_energy[21] # _UCO
        new_z_energy[3:80:4] = z_energy[0:20] + z_energy[20] + z_energy[21] # _UNH_UCO
    else:
        new_z_energy[1:80:4] = z_energy[0:20]
        new_z_energy[2:80:4] = z_energy[0:20]
        new_z_energy[3:80:4] = z_energy[0:20]

    lib = tb.open_file(h5filename, mode='w')
    lib.create_array(lib.root,   'names', new_names)
    lib.create_array(lib.root,'z_energy', new_z_energy)
    lib.root.z_energy._v_attrs.z_min     = z_min
    lib.root.z_energy._v_attrs.z_max     = z_max
    lib.root.z_energy._v_attrs.thickness = thickness
    lib.close()
    return True


def plot_potential_from_h5_file_80type(h5filename):
    lib = tb.open_file(h5filename, mode='r')
    z_energy  = lib.root.z_energy[:]
    z_min     = lib.root.z_energy._v_attrs.z_min
    z_max     = lib.root.z_energy._v_attrs.z_max
    thickness = lib.root.z_energy._v_attrs.thickness
    lib.close()

    # plot
    zpot = np.linspace(z_min, z_max, len(z_energy[0]))
    figure(figsize=(20,10))
    for residx in range(20):
        resname = index_to_restype[residx]

        subplot(4, 5, residx+1)
        axvline(            0, color='k')
        axvline( thickness/2., color='r', linestyle='dashed')
        axvline(-thickness/2., color='r', linestyle='dashed')
        axhline(            0, color='k')

        plot(zpot, z_energy[residx*4  ], color='red'   , label='XXX')
        plot(zpot, z_energy[residx*4+1], color='blue'  , label='XXX_UNH')
        plot(zpot, z_energy[residx*4+2], color='lime'  , label='XXX_UCO')
        plot(zpot, z_energy[residx*4+3], color='purple', label='XXX_UNH_UCO')

        if residx == 0:
            legend(ncol=2,fontsize='large',loc='upper left')

        tick_params(axis='both', which='major', labelsize=15)
        title(resname, fontsize=25)
        ylim( -4, 6)
        xlim(-35,35)
        grid()
    tight_layout()
    savefig('%s.png' % bsnm(h5filename))
    return True


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0] 


def main():
    oldh5 = sys.argv[1]
    newh5 = sys.argv[2]
    names, z_energy, z_min, z_max, thickness = read_potential_h5_file_22type(oldh5)
    write_potential_h5_file_80type(names, z_energy, z_min, z_max, thickness, newh5)
    plot_potential_from_h5_file_80type(newh5)


if __name__ == '__main__':
    main()

