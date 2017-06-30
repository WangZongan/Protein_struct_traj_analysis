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

"s" is the energy_displacement.

Note: energy and s should have the same sign in order to reach minimum at dR = dR0. 


2016.10.21
New Upside changes the write_contact_energies in upside_config.py and struct ContactEnergy in sidechain_radial.cpp.
1. head 'residue1 residue2 r0 width energy' --> 'residue1 residue2 energy distance transition_width', 
   and distance = r0, transition_width = width.
2. CM-CM (center of mass of sidechain) contacts --> CA-CA contacts.

Thus, the default contact type should change to CA. 
'''

__author__  = 'Wang Zongan'
__version__ = '2016-11-18'

import string
import numpy as np
import pandas as pd

import argparse
from itertools import combinations

H_bond=0.88
O_bond=1.24

base_sc_ref = {
     'ALA': np.array([-0.01648328,  1.50453228,  1.20193768]),
     'ARG': np.array([-0.27385093,  3.43874264,  2.24442499]),
     'ASN': np.array([-0.27119135,  2.28878532,  1.32214314]),
     'ASP': np.array([-0.19836569,  2.23864046,  1.36505725]),
     'CYS': np.array([-0.17532601,  1.92513503,  1.34296652]),
     'GLN': np.array([-0.28652696,  2.84800873,  1.60009894]),
     'GLU': np.array([-0.26377398,  2.80887008,  1.69621717]),
     'GLY': np.array([-1.56136239e-02, 5.46052464e-01, -5.67664281e-19]),
     'HIS': np.array([-0.32896151,  2.66635893,  1.42411271]),
     'ILE': np.array([-0.23956042,  2.26489309,  1.49776818]),
     'LEU': np.array([-0.23949426,  2.67123263,  1.3032201 ]),
     'LYS': np.array([-0.26626635,  3.18256448,  1.85836641]),
     'MET': np.array([-0.21000946,  2.79544428,  1.52568726]),
     'PHE': np.array([-0.27214755,  2.83761534,  1.45094383]),
     'PRO': np.array([-1.10993493,  0.89959734,  1.41005877]),
     'SER': np.array([-0.00692474,  1.56683138,  1.475341  ]),
     'THR': np.array([-0.14662723,  1.80061252,  1.42785569]),
     'TRP': np.array([-0.01433503,  3.07506159,  1.56167948]),
     'TYR': np.array([-0.2841611 ,  3.02555746,  1.50123341]),
     'VAL': np.array([-0.02436993,  1.97251406,  1.32782961])}

model_geom = np.zeros((3,3))
model_geom[0] = (-1.19280531, -0.83127186, 0.)  # N
model_geom[1] = ( 0.,          0.,         0.)  # CA
model_geom[2] = ( 1.25222632, -0.87268266, 0.)  # C
model_geom -= model_geom.mean(axis=0)

three_letter_aa = dict(
        A='ALA', C='CYS', D='ASP', E='GLU',
        F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN',
        P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

aa_num = dict([(k,i) for i,k in enumerate(sorted(three_letter_aa.values()))])

one_letter_aa = dict([(v,k) for k,v in three_letter_aa.items()])


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


def calculate_sidechain_CM(N,CA,C,resType):
    '''
    N, CA, C : coordinates of N, CA, C atoms of one given residue
    resType  : 3-letter aa type
    '''
    assert N.shape == CA.shape == C.shape == (1,3)
    z = np.concatenate((N, CA, C),axis=0)
    rot, trans = rmsd_transform(z,model_geom)
    return np.dot(base_sc_ref[resType],rot.T) + trans


def get_coordinates(protein, sel='CA', ref_sidechain_CM=True):
    '''
    protein : prody structure object
    sel     : CA, CB, CM
    '''
    sel_ca    = protein.select('name CA')
    fasta_seq = [s for s in sel_ca.getResnames()]
    nres      = len(fasta_seq)
    
    coords = np.zeros((nres,3))
    hv     = protein.getHierView()
    if sel == 'CA':
        for i,res in enumerate(hv.iterResidues()):
            coords[i] = res.select('name CA').getCoords()
    elif sel == 'CB':
        for i,res in enumerate(hv.iterResidues()):
            at = res.select('name CB') if res.getResname() != 'GLY' else res.select('name CA')
            coords[i] = at.getCoords()
    elif sel == contact == 'CM':
        if ref_sidechain_CM:
            print "Now, use the static reference sidechain CM positions implied by N, CA, and C."
            for i,res in enumerate(hv.iterResidues()):
                coords[i] = calculate_sidechain_CM(res.select('name N').getCoords(),
                                                   res.select('name CA').getCoords(), 
                                                   res.select('name C').getCoords(),
                                                   res.getResname())
        else:
            for i, res in enumerate(hv.iterResidues()):
                if res.getResname() == 'GLY':
                    sel_sd = res.select('name CA')
                else:
                    sel_sd = res.select('sidechain and not name CA') # sidechain includes CB
                coords[i] = prody.calcCenter(sel_sd, weights=sel_sd.getMasses()).round(3)
    else:
        raise ValueError('--contact-type must be either CM or CB or CA')
    return coords


def residue_pair_id_nm_dist(coords, ResName, indices, threshold, distance_cutoff):
    # Get residue index pairs for distance calculation 
        
    res_pairs = np.array([(resi, resj) for (resi, resj) in combinations(indices, 2) if np.absolute(resj-resi) > threshold]) 

    #print indices, res_pairs
    if len(res_pairs) == 0:
        return [], [], []
    else:
        res_pair_distances = compute_distances(coords, res_pairs)
        res_pair_names = np.array([ (ResName[resi], ResName[resj])
            for (resi, resj) in combinations(indices,2) if np.absolute(resj-resi) > threshold])

        res_pair_included          = res_pairs[res_pair_distances < distance_cutoff]
        res_pair_name_included     = res_pair_names[res_pair_distances < distance_cutoff]
        res_pair_distance_included = res_pair_distances[res_pair_distances < distance_cutoff]

        return res_pair_included, res_pair_name_included, res_pair_distance_included


def output_contact_energy_file(f, res_pair_included, res_pair_name_included, res_pair_distance_included, 
                                energy, width, energy_scale, energy_displacement, MJ=None): 
    if MJ is not None: 
        for ip in range(len(res_pair_included)):
            energy_ = MJ[res_pair_name_included[ip][0]][res_pair_name_included[ip][1]] * energy_scale

            print >> f,'%s %s %.4f %.4f %.4f' % (
                        res_pair_included[ip][0], 
                        res_pair_included[ip][1],
                        energy_,
                        res_pair_distance_included[ip] - np.sign(energy_) * energy_displacement,
                        width)
            print >> f,'%s %s %.4f %.4f %.4f' % (
                        res_pair_included[ip][0], 
                        res_pair_included[ip][1],
                        -energy_,
                        res_pair_distance_included[ip] + np.sign(energy_) * energy_displacement,
                        width) 
    else:
        for ip in range(len(res_pair_included)):
            print >> f,'%s %s %.4f %.4f %.4f' % (
                        res_pair_included[ip][0], 
                        res_pair_included[ip][1],
                        energy * energy_scale,
                        res_pair_distance_included[ip] - np.sign(energy) * energy_displacement,
                        width)
            print >> f,'%s %s %.4f %.4f %.4f' % (
                        res_pair_included[ip][0], 
                        res_pair_included[ip][1],
                        -energy * energy_scale,
                        res_pair_distance_included[ip] + np.sign(energy) * energy_displacement,
                        width)


def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
    import re

    if re.match('^([0-9]+(-[0-9]+)?)(,[0-9]+(-[0-9]+)?)*$', s) is None:
        raise argparse.ArgumentTypeError('segments must be of the form 10-30,45,72-76 or similar')

    def parse_seg(x):
        atoms = x.split('-')
        if len(atoms) == 1:
            return np.array([int(atoms[0])])
        elif len(atoms) == 2:
            return np.arange(int(atoms[0]),1+int(atoms[1]))  # inclusive on both ends
        else:
            raise RuntimeError('the impossible happened.  oops.')

    ints = np.concatenate([parse_seg(a) for a in s.split(',')])
    ints = np.array(sorted(set(ints)))   # remove duplicates and sort
    return ints


def compute_distances(coords, pairs): 
    '''
        Input
        ----- ----- 
        coords: (nres, 3), atom coordinates 
        pairs : (N, 2), residue indices 

        --> coords[pairs].shape = (N, 2, 3)
        --> np.diff(coords[pairs],axis=1).shape = (N, 1, 3)
        --> np.sum(np.diff(coords[pairs],axis=1)**2,axis=-1).shape = (N, 1)

        Return
        ----- ----- 
        distances: (N, )
    '''
    return np.sqrt(np.sum(np.diff(coords[pairs],axis=1)**2,axis=-1))[:,0]


def main():
    parser = argparse.ArgumentParser(description='PDB to contact energy file for Upside.',
        usage ="use 'python %(prog)s --help' for more information")
    parser.add_argument('pdb',    help='[required] input .pdb file')
    parser.add_argument('output', help='[required] output file')
    
    parser.add_argument('--width' , default=0.5 , type=float,
        help = 'Parameter of width. Defaulf value is 0.5.')
    parser.add_argument('--energy', default=-3., type=float,
        help = 'Parameter of energy. Defaulf value is -3.')
    parser.add_argument('--energy-scale', default=1., type=float,
        help = 'Parameter of energy scale. Defaulf value is 1.')
    parser.add_argument('--energy-displacement', default=0.5, type=float,
        help = 'Parameter of energy displacement. Defaulf value is 0.5.')

    parser.add_argument('--distance-cutoff' , default=7.5 , type=float,
        help = 'Parameter of distance cutoff: only the distance of a pair of residues that is ' +
               'within this cutoff will be considered as in contact. ' +
               'Defaulf value is 7.5 for CA, which means only the relatively short range contacts will included. ' +
               'Recommended: use 6.5 for CM; 7.0 for CB; 7.5 for CA.')
    parser.add_argument('--contact-type', default='CA', type=str,
        help = 'Three types of contacts are available: CA, CB and CM.\n'+
               'CM refers to the center of mass of sidechain.\n'  +
               'Defaulf value is "CA".')
    parser.add_argument('--primary-sequence-distance-threshold', default=6, type=int,
        help = 'Only when the indices of a pair of residues are larger than the threshold ' +
               'in the primary sequence are those two residues considered in contact. ' +
               'Default value is 6.')

    parser.add_argument('--restraint-group', default=[], action='append', type=parse_segments,
        help = 'If provided, only the residues in the restaint groups will be included in the contact energy calculation. ' +
               'Multiple residue groups may be specified by giving the --restraint-group flag multiple times ' +
               'with different filenames. Note: residue indices starts from 0.' )

    parser.add_argument('--ref-sidechain-CM', default=False, action='store_true',
            help='If turned on, use the static reference sidechain CM positions implied by N, CA, and C.')

    parser.add_argument('--MJ-potential', default=None, type=str,  
        help = 'Path to MJ potential matrix. If turned on, instead of using --energy, MJ potential will used.')

    args = parser.parse_args()

    # atoms selection & obtain coordinates
    import prody
    structure = prody.parsePDB(args.pdb)
    protein   = structure.select('protein')
    sel_ca    = protein.select('name CA')
    ResName   = sel_ca.getResnames()
    fasta_seq = [s for s in ResName]
    nres      = len(fasta_seq)

    # get the coordinate of selected atom of each residue 
    coords = get_coordinates(protein, args.contact_type, args.ref_sidechain_CM)
    
    # calculate contacts & output
    if args.MJ_potential is not None: 
        print "MJ potential is provided. If you don't want to use MJ potential, disable it."
        MJ_df = pd.read_csv(args.MJ_potential, sep='\s+')

    with open(args.output,'w') as f:
        #print >> f, 'residue1 residue2 r0 width energy'
        print >> f, 'residue1 residue2 energy distance transition_width'

        if args.restraint_group:
            print
    #        print 'Restraint groups (uppercase letters are restrained residues)'
            fasta_one_letter = ''.join(one_letter_aa[x] for x in fasta_seq)

            for i, rg in enumerate(args.restraint_group):
                restrained_residues = set(rg)
                assert np.amax(list(restrained_residues)) < len(fasta_seq)
                #print 'group_%i: %s' % (i, ''.join((lett.upper() if i in restrained_residues else lett.lower())
                #                              for i,lett in enumerate(fasta_one_letter)))

                rp_included, rp_nm_included, rp_dist_included = residue_pair_id_nm_dist(coords, ResName, list(restrained_residues), 
                            args.primary_sequence_distance_threshold, args.distance_cutoff)

                output_contact_energy_file(f, rp_included, rp_nm_included, rp_dist_included, 
                                        args.energy, args.width, args.energy_scale, args.energy_displacement, MJ_df)
        else:
            rp_included, rp_nm_included, rp_dist_included = residue_pair_id_nm_dist(coords, ResName, range(nres), 
                    args.primary_sequence_distance_threshold, args.distance_cutoff)

            output_contact_energy_file(f, rp_included, rp_nm_included, rp_dist_included,
                                    args.energy, args.width, args.energy_scale, args.energy_displacement, MJ_df)


if __name__ == '__main__':
    main()
