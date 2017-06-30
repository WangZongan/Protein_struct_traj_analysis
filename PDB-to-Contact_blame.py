#!/usr/bin/env python
'''
Generate a text file that contains information for a contact energy function for UPSIDE.

The first line of the file should be a header: "residue1 residue2 r0 width energy".
The remaining lines should contain space separated values.
Example file:
residue1 residue2 r0 width energy
0       1       5.9336  2.0000  -3.0000
0       3       4.5248  2.0000  -3.0000
0       146     5.6584  2.0000  -3.0000
....

The formula of the contact energy:

                                energy
    contact_energy = ---------------------------- , in which
                                1
                     1 + exp[ ----- * (dR - dR0)]
                              width

    dR  = |R(resi) - R(resj)|, and
    dR0 = |R_target(resi) - R_target(resj)|.

R(resi) is the center of mass of residue i's sidechain or CB of residue i.


2016.07.27
Using two entries for each pair of contact --> contact energy reaches minimum when at dR = dR0. 

                                    energy                             - energy
==> contact_energy = ---------------------------------- + ---------------------------------- 
                                1                                    1
                     1 + exp[ ----- * (dR - (dR0 - s))]   1 + exp[ ----- * (dR - (dR0 + s))]
                              width                                width

Note: energy and s should have the same sign in order to reach minimum at dR = dR0. 


2016.10.21
New Upside changes the write_contact_energies in upside_config.py and struct ContactEnergy in sidechain_radial.cpp.
1. head 'residue1 residue2 r0 width energy' --> 'residue1 residue2 energy distance transition_width', 
   and distance = r0, transition_width = width.
2. CM-CM (center of mass of sidechain) contacts --> CA-CA contacts.

Thus, the default contact type should change to CA. 
'''

__author__  = 'Wang Zongan'
__version__ = '2016-12-15'

import os
import string
import numpy as np
import pandas as pd

import argparse
from itertools import combinations

import mdtraj as md

H_bond=0.88
O_bond=1.24

three_letter_aa = dict(
        A='ALA', C='CYS', D='ASP', E='GLU',
        F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN',
        P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

aa_num = dict([(k,i) for i,k in enumerate(sorted(three_letter_aa.values()))])

one_letter_aa = dict([(v,k) for k,v in three_letter_aa.items()])


def read_contact_energy_file(fpath):
    print 'Load contact energy file.'
    df = pd.read_csv(fpath, delim_whitespace=True)

    pairs = np.column_stack((df.residue1.as_matrix()[::2], df.residue2.as_matrix()[::2])).astype(np.int16)
    energy = df.energy.as_matrix().astype(np.float32)
    distance = df.distance.as_matrix().astype(np.float32)
    transition_width = df.transition_width.as_matrix().astype(np.float32)

    return pairs, energy, distance, transition_width


def load_mdtraj(pdbpath):
    print 'Load PDB.'
    return md.load(pdbpath)


def traj_compute_contacts(traj, pairs):
    return md.compute_contacts(traj, pairs, scheme='ca')[0]*10 # unit: A


def compute_contact_energy(distances, ref_distance, energy, transition_width):
    '''
    Input
    ----------
    distances        : (nframes, npairs)
    ref_distance     : (2 * npairs, ) 
    energy           : (2 * npairs, )
    transition_width : (2 * npairs, )

    Output
    ----------
    energies         : (nframes, )

    
    For each term:
                                energy         
    contact_energy = ---------------------------- 
                                1                              
                     1 + exp[ ----- * (dR - dR0)]
                              width             
    '''
    print 'Blame contact energy at each step.'

    compute_ = lambda e, w, dr, dr0: np.sum(e/(1+np.exp((dr-dr0)/w)), axis=1)

    return compute_(energy[0::2], transition_width[0::2], distances, ref_distance[0::2]) + \
           compute_(energy[1::2], transition_width[1::2], distances, ref_distance[1::2])
    

bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def main():
    parser = argparse.ArgumentParser(description='PDB to contact energy file for Upside.',
        usage ="use 'python %(prog)s --help' for more information")
    parser.add_argument('pdb',    help='[required] input .pdb file')
    parser.add_argument('energy', help='[required] input contact energy file')
    args = parser.parse_args()

    pairs, energy, distance, transition_width = read_contact_energy_file(args.energy)
    traj = load_mdtraj(args.pdb)

    distances = traj_compute_contacts(traj, pairs) 
    energies  = compute_contact_energy(distances, distance, energy, transition_width)

    with open('%s.%s.contact_energy.dat' % (bsnm(args.pdb), bsnm(args.energy)), 'w') as f:
        for e in energies:
            print >> f, '%.3f' % e
    
    
def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


if __name__ == '__main__':
    import time
    sta = time.time()
    main()
    print '\nrunning time: %s' % sec_to_hr_min_sec(time.time() - sta)


