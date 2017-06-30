#!/usr/bin/env python
'''
Extract output array from upside-generated h5 files.
'''

__author__  = 'Wang Zongan'
__version__ = '2016-09-13'

import os
import sys
import string
import numpy as np
import tables as tb


def read_traj(s, path):
    '''
    From John Jumper's energy_blame.py
    '''
    d=dict()
    with tb.open_file(path) as t:
        o = t.root.output
        d['seq']    = t.root.input.sequence[:]
        d['seq'][d['seq']=='CPR'] = 'PRO'
        print 'n_res', len(d['seq'])
        d['pos']    = o.pos[s:,0]
        d['pot']    = o.potential[s:,0]
        d['strain'] = o.rotamer_1body_energy0[s:]
        d['cov']    = o.rotamer_1body_energy1[s:]
        d['hydro']  = o.rotamer_1body_energy2[s:]
        d['sc']     = o.rotamer_free_energy  [s:]
        d['pair']   = d['sc'] - d['strain'] - d['cov'] #- d['env']
        d['rama']   = o.rama_map_potential[s:]
        d['hb']     = read_hb(t)[s:]
        d['Rg']     = np.sqrt(np.var(d['pos'],axis=1).sum(axis=-1))
        phe = t.root.input.potential.hbond_energy._v_attrs.protein_hbond_energy
        d['hb_energy'] = phe*d['hb'][...,0].sum(axis=-1)
        d['env']    = o.nonlinear_coupling[s:]
    return d


def read_hb(tr):
    '''
    From John Jumper's energy_blame.py
    '''
    n_res = tr.root.input.pos.shape[0]/3
    don_res =  tr.root.input.potential.infer_H_O.donors.id[:,1] / 3
    acc_res = (tr.root.input.potential.infer_H_O.acceptors.id[:,1]-2) / 3

    n_hb = tr.root.output.hbond.shape[1]
    hb_raw   = tr.root.output.hbond[:]
    hb = np.zeros((hb_raw.shape[0],n_res,2,2))
    
    hb[:,don_res,0,0] =    hb_raw[:,:len(don_res)]
    hb[:,don_res,0,1] = 1.-hb_raw[:,:len(don_res)]

    hb[:,acc_res,1,0] =    hb_raw[:,len(don_res):]
    hb[:,acc_res,1,1] = 1.-hb_raw[:,len(don_res):]

    return hb



def main():
    h5_file_path = sys.argv[1]
    stride = int(sys.argv[2])

    #output_array = lib.get_node('/output', array_name)
    lib = tb.open_file(h5_file_path)
    pos = lib.root.output.pos[:,0] 
    pot = lib.root.output.potential[:,0]
    hbond = np.around(np.sum(lib.root.output.hbond[:], axis=1),1)

    Rg = np.sqrt(np.var(pos,axis=1).sum(axis=-1))

    output_file_name = os.path.splitext(os.path.basename(h5_file_path))[0]
    with open(output_file_name+'.potential.dat','w') as f:
        for val in pot[::stride]:
            print >> f, val 
    with open(output_file_name+'.hbond.dat','w') as f:
        for val in hbond[::stride]:
            print >> f, val
    with open(output_file_name+'.Rg.dat','w') as f:
        for val in Rg[::stride]:
            print >> f, val

    lib.close()

if __name__ == '__main__':
    main()


