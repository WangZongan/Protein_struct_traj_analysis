#!/usr/bin/python
'''
Given a pdb file, up to user's choice, extract the residue properties of the protein.
'''

__author__  = 'Wang Zongan'

import os
import string, math
import numpy as np
import collections
import tables as tb

import mdtraj as md

oneletter_threeletter = dict(
        A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

threeletteraa2number = dict([(oneletter_threeletter[k],i) for i,k in enumerate(sorted(oneletter_threeletter.keys()))])
number2threeletteraa = dict([(i,oneletter_threeletter[k]) for i,k in enumerate(sorted(oneletter_threeletter.keys()))])

backbone   = ['C','CA','N']
backbone_O = ['C','CA','N','O']


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def calc_sasa_sidechain_nh_co(mdtraj_protein):
    top = mdtraj_protein.topology
    sasa_atom = md.shrake_rupley(mdtraj_protein, probe_radius=0.14, n_sphere_points=960, mode='atom')[0]
    assert top.n_atoms == len(sasa_atom)

    res_sasa_sc = np.zeros(top.n_residues)
    res_sasa_nh = np.zeros(top.n_residues)
    res_sasa_co = np.zeros(top.n_residues)

    for ires in range(top.n_residues):
        res = top.residue(ires)
        atnm = [res.atom(i).name for i in range(res.n_atoms)]

        is_cter = 'OXT' in atnm
        is_nter = 'H2' in atnm

        # sc
        if res.name != 'GLY':
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.is_sidechain:
                    res_sasa_sc[ires] += sasa_atom[at.index]
        else:
            h2_sasa = 0; h3_sasa = 0
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'H2':
                    h2_sasa = sasa_atom[at.index]
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'H3':
                    h3_sasa = sasa_atom[at.index]
            res_sasa_sc[ires] = (h2_sasa+h3_sasa)/2.

        # co
        if is_cter:
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'C': c_sasa = sasa_atom[at.index]
                if at.name == 'O': o_sasa = sasa_atom[at.index]
                if at.name == 'OXT': oxt_sasa = sasa_atom[at.index]
            res_sasa_co[ires] = o_sasa+oxt_sasa+c_sasa
        else:
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'C': c_sasa = sasa_atom[at.index]
                if at.name == 'O': o_sasa = sasa_atom[at.index]
            res_sasa_co[ires] = o_sasa+c_sasa

        # nh
        if is_nter:
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'N': res_sasa_nh[ires]   += sasa_atom[at.index]
                if at.name == 'H': res_sasa_nh[ires]   += sasa_atom[at.index]
                if at.name == 'H2': res_sasa_nh[ires]  += sasa_atom[at.index]
                if at.name == 'H3': res_sasa_nh[ires]  += sasa_atom[at.index]
        else:
            for nat in range(res.n_atoms):
                at = res.atom(nat)
                if at.name == 'N': n_sasa   = sasa_atom[at.index]
                if at.name == 'H':
                    h_sasa = sasa_atom[at.index]
                else:
                    h_sasa = 0
            res_sasa_nh[ires] = h_sasa+n_sasa

    return res_sasa_sc, res_sasa_nh, res_sasa_co


#### <----- ----- ----- ----- ----- ----- ----- -----> ####
output_head_string  = '#Index restype secseq '
output_head_string += 'sasa_sc sasa_nh sasa_co '
output_head_string += 'CBx CBy CBz'


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Prepare input file',
        usage ='use "%(prog)s --help" for more information')

    parser.add_argument('pdb',    help='[required] input .pdb file')
    parser.add_argument('--membrane-thickness', default=30.0, type=float,
	help='Thickness of membrane bilayer, if the protein is a membrane protein.')

    args = parser.parse_args()

    protein = md.load_pdb(args.pdb)
    top     = protein.topology
    nres    = top.n_residues

    restypes  = [top.residue(i).name for i in range(nres)]

    '''
    = = = = = = = = = = = = = = = = = = = =
    coordinates
    = = = = = = = = = = = = = = = = = = = =
    '''
    CB_coords = np.zeros((nres,3))
    for i in range(nres):
        res = top.residue(i)
        if res.name != 'GLY':
            for j in range(res.n_atoms):
                at = res.atom(j)
                if at.name == 'CB':
                    CB_coords[i] = protein.xyz[0][at.index]
        else:
            for j in range(res.n_atoms):
                at = res.atom(j)
                if at.name == 'CA':
                    CB_coords[i] = protein.xyz[0][at.index]


    '''
    = = = = = = = = = = = = = = = = = = = =
    DSSP
    = = = = = = = = = = = = = = = = = = = =
    Including Secondary Structure, SASA, and relative accessibility.
    Secondary structure given by DSSP definition.

    DSSP codes in Bio.PDB
    DSSP Code   Secondary structure
    H   alpha-helix
    B   Isolated beta-bridge residue
    E   beta-Strand
    G   3-10 helix
    I   Pi-helix
    T   Turn
    S   Bend (region of high curvature)
    -   Other
    C   coil
    N   non-structure.

    Solvent assessible surface area, calculated by DSSP algorithm.
    Relative accessiblity, calculated by DSSP algorithm.
    = = = = = = = = = = = = = = = = = = = =
    '''
    secseq = md.compute_dssp(protein, simplified=False)[0]
    secseq = np.where(secseq == ' ', 'N', secseq)

    res_sasa_sc, res_sasa_nh, res_sasa_co = calc_sasa_sidechain_nh_co(protein)

    '''
    = = = = = = = = = = = = = = = = = = = =
    output
    = = = = = = = = = = = = = = = = = = = =
    '''
    with open(bsnm(args.pdb) + '.residue_CB_SASA.dat','w') as f:
        print >> f, '#thicnkess %.1f' % args.membrane_thickness
        print >> f, output_head_string
        for resnum in range(nres):
            print >> f, '%4i %3s %1s %6.3f %6.3f %6.3f %8.3f %8.3f %8.3f' % \
                (resnum, restypes[resnum], secseq[resnum],
                res_sasa_sc[resnum], res_sasa_nh[resnum], res_sasa_co[resnum],
                CB_coords[resnum][0], CB_coords[resnum][1], CB_coords[resnum][2])


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    end_time   = time.time()
    print
    print 'Total running time:'
    print "--- %s seconds ---" % round((end_time - start_time),3)
