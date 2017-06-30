#!/usr/bin/env python
'''
This program is used for calculating and plotting the contact map for a given pdb file.

The distance between two residues can be defined in different ways:
1. the distance between a specific pair of atoms (i.e., CA-CA or CB-CB),
2. the shortest distance among the atoms belonging to residue i and those belonging to residue j
3. the distance between the centers of mass of the two residues
'''

__author__  = 'Wang Zongan'
__version__ = '2016-05-02'

import os
import sys 
import string
import numpy as np
import cPickle as cp

import Bio.PDB

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *

backbone = ['C','CA','N','O']

def residue_distance_sidechainCM(res_1,res_2):
    def residue_sidechainCM(res):
        atoms = []
        for atom in res.get_atom():
            if atom not in backbone:
                atoms.append(atom)
        coord = np.zeros(3)
        mass = 0
        for atom in atoms:
            coord += atom.get_coord()*atom.mass
            mass += atom.mass
        return coord/mass
    res_1_coord = residue_sidechainCM(res_1)
    res_2_coord = residue_sidechainCM(res_2)
    diff_vec = res_1_coord - res_2_coord
    return np.sqrt(np.sum(diff_vec*diff_vec,axis=0))


def find_contact_distance_cutoff(structure, contact_type='CA', contact_cutoff=7.5):
    nmodel = len(structure)
    # exclude het-residues and water molecules
    nres = len([res for res in structure[0].get_residues() if res.get_id()[0] == ' '])

    '''
    <Yuan et al. BMC Bioinformatics 13 (2012) 292>
    contact_type:
    1. 'CA': CA-CA contact, default cutoff 7.5
    2. 'CB': CB-CB contact, default cutoff 7.0
    3. 'CM': sidechain CM-CM contact, default cutoff 6.5
    '''
    contacts = np.zeros((nmodel,nres,nres))
    dists    = np.zeros((nmodel,nres,nres)) 
    for nm, model in enumerate(structure):
        residues = [res for res in structure[nm].get_residues() if res.get_id()[0] == ' '] 
        assert len(residues) == nres
        for r1 in range(nres-1):
            res_1 = residues[r1]
            for r2 in range(r1+1,nres):
                res_2 = residues[r2]
                if (r2 - r1) >= 4:
                    if contact_type == 'CA':
                        dist = res_1['CA'] - res_2['CA']
                    elif contact_type == 'CB':
                        res_1_atom = res_1['CB'] if res_1.get_resname() != 'GLY' else res_1['CA']
                        res_2_atom = res_2['CB'] if res_2.get_resname() != 'GLY' else res_2['CA']
                        dist = res_1_atom - res_2_atom
                    elif contact_type == 'CM':
                        dist = residue_distance_sidechainCM(res_1,res_2)
                    if dist <= contact_cutoff:
                        contacts[nm,r1,r2] = 1
                        contacts[nm,r2,r1] = 1
                        dists[nm,r1,r2]    = dist 
    return contacts, dists 


def plot_contact_map(contacts, nres, pdbfile, contact_type, cutoff, spacing=10, 
                plot_color_map='Reds', plot_format='png', ref_contacts=None):  

    def make_ticklabels(nres,spacing):
        xtl = []; ytl = []
        for i in range(nres):
            if i % spacing == 0:
                xtl.append(str(i)); ytl.append(str(i))
            else:
                xtl.append(''); ytl.append('')
        return xtl,ytl
    xtl,ytl = make_ticklabels(nres,spacing)

    Z = np.sum(contacts,axis=0)/contacts.shape[0] # average contact map of the trajectory

    if ref_contacts is not None:
        ref_Z = np.sum(ref_contacts,axis=0)/ref_contacts.shape[0]
        assert ref_Z.shape == Z.shape
        # set the lower half of Z to be the lower half of ref_Z
        row, col = Z.shape 
        for irow in range(row): 
            for icol in range(row)[0:irow]: 
                Z[irow,icol] = ref_Z[irow,icol]   

    fig = figure(figsize=(12,10))
    suptitle('%s Contact Map: cutoff=%s'%(contact_type, str(cutoff)),fontsize=25,x=0.435,y=0.98)
    title(pdbfile, fontsize=20)
    my_cmap = get_cmap(plot_color_map)
    my_cmap.set_under('w')

    sns.heatmap(Z, vmin=0.00001, vmax=1, cmap=my_cmap, linewidths=0.001, linecolor='white',
              xticklabels=xtl, yticklabels=ytl)

    tick_params(axis='both', which='major', labelsize=15)
    fig.text(0.4, 0.03, 'Residue', ha='center', va='center', fontsize=25)
    fig.text(0.02, 0.5, 'Residue', ha='center', va='center', fontsize=25, rotation='vertical')
    grid()
    tight_layout(rect=[0.03, 0.05, 1, 0.95]) # default is (0, 0, 1, 1) [left, bottom, right, top]
    savefig('%s.contact_map.%s.%s.%s'%(pdbfile, contact_type, str(cutoff), plot_format), format=plot_format)  


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Calculate and plot contact map for a given pdb file, ' +
                      'which can contain a trajectory of structures.',
        usage ='use "python %(prog)s --help" for more information')
    parser.add_argument('pdb', help='[required] input PDB file')
    parser.add_argument('--contact-type', default='CA', type=str,
        help = 'Contact type, options are CA, CB, CM, which stand for CA-CA contact, '+
               'CB-CB contact, and sidechain CM-CM contact, respectively. '+
               'Now all the contact types supported are distance cutoff contacts.' +
               'The default value is CA.')
    parser.add_argument('--contact-cutoff', default=7.5, type=float,
        help = 'Suggested values for CA, CB, and CM contacts are 7.5, 7.0, 6.5 A, respectively. ' +
               'The default value is 7.5 for the default CA contact.')

    parser.add_argument('--reference-model', default=None, type=str,
        help = 'If provided, the output plot of contact map will contain 2 parts, ' + 
               'the lower half being the contact map of the reference and ' + 
               'the upper half being the contact map of the pdb structure.')

    parser.add_argument('--plot-tick-spacing', default=10, type=int,
        help = 'Used for adjust the tick labels on the output contact map. The default value is 10.')
    parser.add_argument('--plot-color-map', default='Reds', type=str, 
        help = "Color map used in plotting, spectral by default. " +
               "Any color map supported by matplotlib can be used. " + 
               "Examples are: 'Blues', 'GnBu', 'BrBG', 'gist_rainbow', etc. " + 
               "(Ref: http://matplotlib.org/xkcd/examples/color/colormaps_reference.html)")
    parser.add_argument('--plot-format', default='png', type=str,
        help = 'Format of output plot, PNG format by default. ' +
               'Any format supported by matplotlib can be used.')

    parser.add_argument('--output-contact-energy-file', default=None, type=str,
        help = 'If turned on, output the file for contact energy for upside.')
    parser.add_argument('--output-contact-energy-file-freq-cutoff', default=0.5, type=float,
        help = 'Output contact pairs with frequency larger than the cutoff. ' +
               'The default value is 0.5.')
    args = parser.parse_args()

    pdbfile = args.pdb
    struct = Bio.PDB.PDBParser(QUIET=True).get_structure('protein',pdbfile)
    nres = len([res for res in struct[0].get_residues() if res.get_id()[0] == ' '])
    nmodel = len(struct)

    print '=' * 40
    print '%i models detected in the pdb file.' % nmodel
    print '%i residues in the model' % nres
    print '=' * 40

    contacts, dists = find_contact_distance_cutoff(struct, args.contact_type, args.contact_cutoff)

    if args.reference_model is not None:
        ref_structure = Bio.PDB.PDBParser(QUIET=True).get_structure('protein',args.reference_model)
        ref_contacts, _ = find_contact_distance_cutoff(ref_structure, args.contact_type, args.contact_cutoff)
        plot_contact_map(contacts, nres, os.path.basename(pdbfile), 
            contact_type=args.contact_type, cutoff=args.contact_cutoff, spacing=args.plot_tick_spacing, 
            plot_color_map=args.plot_color_map, plot_format=args.plot_format, ref_contacts=ref_contacts) 
    else:
        plot_contact_map(contacts, nres, os.path.basename(pdbfile), 
            contact_type=args.contact_type, cutoff=args.contact_cutoff, spacing=args.plot_tick_spacing, 
            plot_color_map=args.plot_color_map, plot_format=args.plot_format) 

    if args.output_contact_energy_file is not None:
        '''
        The formula of the contact energy: 

                                  energy 
        contact_energy = ---------------------------- , in which 
                                    1                     
                         1 + exp[ ----- * (dR - dR0)] 
                                  width 

        dR  = |R(resi) - R(resj)|, and 
        dR0 = |R_target(resi) - R_target(resj)|. 

        Example file: 
        residue1 residue2 r0 width energy
        0       1       5.9336  2.0000  -3.0000
        0       3       4.5248  2.0000  -3.0000
        0       146     5.6584  2.0000  -3.0000
        ....

        r0's are set to be the average distance of the residue pairs in contact.
        Energies are set to be -3.0*contact_frequency. 
        '''
        average_dists = np.sum(dists,axis=0)/nmodel
        average_contacts = np.sum(contacts,axis=0)/nmodel
        f = open(args.output_contact_energy_file,'w')
        print >> f,'residue1 residue2 r0 width energy'
        for i in range(nres):
            for j in range(nres):
                if average_dists[i,j] > 0:
                    if average_contacts[i,j] >= args.output_contact_energy_file_freq_cutoff: 
                        print >> f,'%i\t%i\t%.4f\t%.4f\t%.4f' % (i,j,average_dists[i,j],
                                                              2.0/average_contacts[i,j],
                                                             -3.0*average_contacts[i,j])
        f.close()
        
if __name__ == '__main__':
    main()
