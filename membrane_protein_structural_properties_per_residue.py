#!/usr/bin/python
'''
Given a pdb file, up to user's choice, extract the residue properties of the protein.
'''

__author__  = 'Wang Zongan'
__version__ = '2017.02.17'

import os
import string, math
import numpy as np
import collections 
import simplejson as json
import jsonpickle as jp
import tables as tb

import Bio.PDB  # hierarchy: structure --> model --> chain --> residue --> atom 
from Bio.PDB.Vector import calc_angle 

oneletter_threeletter = dict(
        A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN', P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

threeletteraa2number = dict([(oneletter_threeletter[k],i) for i,k in enumerate(sorted(oneletter_threeletter.keys()))])
number2threeletteraa = dict([(i,oneletter_threeletter[k]) for i,k in enumerate(sorted(oneletter_threeletter.keys()))])

backbone   = ['C','CA','N'] 
backbone_O = ['C','CA','N','O'] 

nres_in_model = lambda model: len([res for res in model.get_residues() if res.get_id()[0] == ' ']) 

bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def get_sidechain_com(residue):

    sidechain_atoms = []
    for atom in residue.get_atom():
        if atom.get_name() in backbone_O: continue 
        elif atom.get_name() in ['HA', 'HA1', 'H']: continue 
        else:
            sidechain_atoms.append(atom) 
    if sidechain_atoms == []: 
        print 'No atoms in sidechain detected !'

    coord = np.zeros(3)
    mass  = 0 
    for atom in sidechain_atoms:
        coord += atom.get_coord()*atom.mass
        mass += atom.mass

    assert mass != 0
    return coord/mass


#### <----- ----- ----- ----- hydrogen bond ----- ----- ----- -----> ####
'''
    Geometric criteria of hbond identification given by papers:
    <Protein Eng. 10 (1997) 999>, <Protein Eng. 15 (2002) 359>

    D2-D1-D-H...A-A1
    
    Rule 1: D-A     <= 3.9 A;
    Rule 2: H-A     <= 2.5 A;
    Rule 3: D-H-A   > 90.0 degree;
    Rule 4: A1-A-H  > 90.0 degree;
    Rule 5: A1-A-D  > 90.0 degree.

    Geometric criteria of ideal Calpha-H hbond:
    <PNAS 98(2001)9056>

    D2-D1-CA-H...O-A1

    1. Calpha-O       <= 3.8 A;
    2. H-O            <= 2.7 A; 
    3. Calpha-H-O      = 180 degree;
    4. A1-O-H          = 120 degree
    5. elevation angle =   0 degreen (vector<Calpha-H> & amide plane) 
    
    Bifurcated hbonds exist. 

    H-bond's types: bb-bb, bb-sc, sc-sc, HA-bb, HA-sc 

    aa   donors                                 acceptor 
    ALA  HN, H(HA)                              O 
    ARG  HN, H(HA), HE, 1HH1, 2HH1, 1HH2, 2HH2  O
    ASN  HN, H(HA), 1HD2, 2HD2                  O, OD1
    ASP  HN, H(HA)                              O, OD1, OD2
    CYS  HN, H(HA)                              O
    GLN  HN, H(HA), 1HE2, 2HE2                  O, OE1
    GLU  HN, H(HA)                              O, OE1, OE2
    GLY  HN, H(HA1), HA2                        O 
    HIS  HN, H(HA), HE2                         O 
    ILE  HN, H(HA)                              O
    LEU  HN, H(HA)                              O 
    LYS  HN, H(HA), 1HZ, 2HZ, 3HZ               O 
    MET  HN, H(HA)                              O 
    PHE  HN, H(HA)                              O 
    PRO      H(HA)                              O 
    SER  HN, H(HA), HG                          O, OG 
    THR  HN, H(HA), HG1                         O, OG1 
    TRP  HN, H(HA), 1HB, 2HB, HE1               O 
    TYR  HN, H(HA), HH                          O, OH  
    VAL  HN, H(HA)                              O 
'''


normalize = lambda vec: vec/np.sqrt(np.sum(vec*vec,axis=-1))
calc_vec_length = lambda vec: np.sqrt(np.sum(vec*vec,axis=-1))


def indeed_hbonded(D_atom,H_atom,A_atom,A1_atom):  
    '''
    D-H...A-A1
    
    Rule 1: D-A     <= 3.9 A;
    Rule 2: H-A     <= 2.7 A;
    Rule 3: D-H-A   > 90.0 degree;
    Rule 4: A1-A-H  > 90.0 degree;
    Rule 5: A1-A-D  > 90.0 degree.
    '''    
    if D_atom - A_atom <= 3.9:
        if H_atom - A_atom <= 2.7:
            if np.rad2deg(calc_angle(D_atom.get_vector(),H_atom.get_vector(),A_atom.get_vector())) > 90:
                if np.rad2deg(calc_angle(A1_atom.get_vector(),A_atom.get_vector(),H_atom.get_vector())) > 90: 
                    if np.rad2deg(calc_angle(A1_atom.get_vector(),A_atom.get_vector(),D_atom.get_vector())) > 90: 
                        return True 
                    else: 
                        return False 
                else:
                    return False 
            else:
                return False 
        else:
            return False 
    else:
        return False 


hbond_donor_related_atoms_dict = { 
    'ALA':[['HN','N'],['H','CA']],               
    'ARG':[['HN','N'],['H','CA'],['HE','NE'],['1HH1','NH1'],['2HH1','NH1'],['1HH2','NH2'],['2HH2','NH2']], 
    'ASN':[['HN','N'],['H','CA'],['1HD2','ND2'],['2HD2','ND2']], 
    'ASP':[['HN','N'],['H','CA']], 
    'CYS':[['HN','N'],['H','CA']],               
    'GLN':[['HN','N'],['H','CA'],['1HE2','NE2'],['2HE2','NE2']], 
    'GLU':[['HN','N'],['H','CA']],               
    'GLY':[['HN','N'],['H','CA'],['HA2','CA']], 
    'HIS':[['HN','N'],['H','CA'],['HE2','NE2']],         
    'ILE':[['HN','N'],['H','CA']], 
    'LEU':[['HN','N'],['H','CA']],               
    'LYS':[['HN','N'],['H','CA'],['1HZ','NZ'],['2HZ','NZ'],['3HZ','NZ']], 
    'MET':[['HN','N'],['H','CA']],               
    'PHE':[['HN','N'],['H','CA']], 
    'PRO':[           ['H','CA']],               
    'SER':[['HN','N'],['H','CA'],['HG','OG']], 
    'THR':[['HN','N'],['H','CA'],['HG1','OG1']],         
    'TRP':[['HN','N'],['H','CA'],['1HB','CB'],['2HB','CB'],['HE1','NE1']], 
    'TYR':[['HN','N'],['H','CA'],['HH','OH']],          
    'VAL':[['HN','N'],['H','CA']]
}

hbond_acceptor_related_atoms_dict = {
    'ALA':[['O','C','CA']],       
    'ARG':[['O','C','CA']],     
    'ASN':[['O','C','CA'], ['OD1','CG','CB']],       
    'ASP':[['O','C','CA'], ['OD1','CG','CB'], ['OD2','CG','CB']], 
    'CYS':[['O','C','CA']],       
    'GLN':[['O','C','CA'], ['OE1','CD','CG']], 
    'GLU':[['O','C','CA'], ['OE1','CD','CG'], ['OE2','CD','CG']], 
    'GLY':[['O','C','CA']],
    'HIS':[['O','C','CA']],    
    'ILE':[['O','C','CA']],     
    'LEU':[['O','C','CA']],          
    'LYS':[['O','C','CA']],
    'MET':[['O','C','CA']],    
    'PHE':[['O','C','CA']],    
    'PRO':[['O','C','CA']],          
    'SER':[['O','C','CA'], ['OG','CB','HG']], 
    'THR':[['O','C','CA'], ['OG1','CB','HG1']], 
    'TRP':[['O','C','CA']],       
    'TYR':[['O','C','CA'], ['OH','CZ','HH']], 
    'VAL':[['O','C','CA']]
}


def find_hbond_pairs_in_two_residues(residue_1,residue_2): 
    hbond_pairs = [] # every element is [donor_resid, acceptor_resid, donor_related_atoms, acceptor_related_atoms]
    # residue 1 --> acceptor 
    for acceptor_related_atoms in hbond_acceptor_related_atoms_dict[residue_1.get_resname()]:
        for donor_related_atoms in hbond_donor_related_atoms_dict[residue_2.get_resname()]:
            try:
                D_atom  = residue_2[donor_related_atoms[1]]
                H_atom  = residue_2[donor_related_atoms[0]]
                A_atom  = residue_1[acceptor_related_atoms[0]]
                A1_atom = residue_1[acceptor_related_atoms[1]]
                if indeed_hbonded(D_atom, H_atom, A_atom, A1_atom):
                    hbond_pairs.append([residue_2.get_full_id()[3], residue_1.get_full_id()[3], 
                                        donor_related_atoms, acceptor_related_atoms])
            except KeyError:
                continue 
    # residue 2 --> acceptor 
    for acceptor_related_atoms in hbond_acceptor_related_atoms_dict[residue_2.get_resname()]:
        for donor_related_atoms in hbond_donor_related_atoms_dict[residue_1.get_resname()]:
            try:
                D_atom  = residue_1[donor_related_atoms[1]]
                H_atom  = residue_1[donor_related_atoms[0]]
                A_atom  = residue_2[acceptor_related_atoms[0]]
                A1_atom = residue_2[acceptor_related_atoms[1]]
                if indeed_hbonded(D_atom, H_atom, A_atom, A1_atom): 
                    hbond_pairs.append([residue_1.get_full_id()[3], residue_2.get_full_id()[3], 
                                        donor_related_atoms, acceptor_related_atoms])
            except KeyError:
                continue 
    return hbond_pairs


donor_info_tuple = collections.namedtuple('donor',
    'chain_id res_index_in_chain res_index_in_model resname donor_atom_name'.split())
donor_related_atom_names_tuple  = collections.namedtuple('donor_related_atom_names', 'D H'.split()) 
donor_related_atom_coords_tuple = collections.namedtuple('donor_related_atom_coords', 'D H'.split()) 

acceptor_info_tuple = collections.namedtuple('acceptor',
    'chain_id res_index_in_chain res_index_in_model resname acceptor_atom_name'.split())  
acceptor_related_atom_names_tuple  = collections.namedtuple('acceptor_related_atom_names', 'A A1 A2'.split())  
# A2 is connected to either A1 or A 
acceptor_related_atom_coords_tuple = collections.namedtuple('acceptor_related_atom_coords', 'A A1 A2'.split()) 


donor_info_string    = 'donor_chain_id donor_res_index_in_chain donor_res_index_in_model donor_resname '
acceptor_info_string = 'acceptor_chain_id acceptor_res_index_in_chain acceptor_res_index_in_model acceptor_resname '
hbond_tuple_string   = donor_info_string + 'D_name H_name D_coord H_coord '
hbond_tuple_string  += acceptor_info_string + 'A_name A1_name A2_name A_coord A1_coord A2_coord'
hbond_tuple          = collections.namedtuple('HBond', hbond_tuple_string.split()) 


def find_hbond_in_model(model): 
    hbonds_pairs_in_model = []

    nres = nres_in_model(model)

    # simple counts just to see to be or not to be 
    simple_NH_hbond_count = np.zeros(nres,dtype=int)  
    simple_CO_hbond_count = np.zeros(nres,dtype=int) 
    simple_SC_hbond_count = np.zeros(nres,dtype=int) 
    simple_CA_hbond_count = np.zeros(nres,dtype=int)

    residues_in_model = [residue for residue in model.get_residues()]

    for resnum_1 in range(nres-1):
        residue_1 = residues_in_model[resnum_1]
        resname_1 = residue_1.get_resname() 
        resid_1   = residue_1.get_full_id()[3]  # (" ", 10, "A") 
        
        for resnum_2 in range(resnum_1+1,nres):
            residue_2 = residues_in_model[resnum_2]
            resname_2 = residue_2.get_resname()
            resid_2   = residue_2.get_full_id()[3]

            hbond_pairs = find_hbond_pairs_in_two_residues(residue_1, residue_2)
            if hbond_pairs != []:  # [donor_resid, acceptor_resid, donor_related_atoms, acceptor_related_atoms]
                for hbond_pair in hbond_pairs:
                    try:
                        if hbond_pair[0] == resid_1: # res 1 is the donor                         
                            if hbond_pair[2][0] == 'HN': 
                                simple_NH_hbond_count[resnum_1] += 1 
                            elif hbond_pair[2][0] == 'H': 
                                simple_CA_hbond_count[resnum_1] += 1
                            else: 
                                simple_SC_hbond_count[resnum_1] += 1

                            if hbond_pair[3][0] == 'O': 
                                simple_CO_hbond_count[resnum_2] += 1
                            else: 
                                simple_SC_hbond_count[resnum_2] += 1

                            hbonds_pairs_in_model.append(
                                hbond_tuple(
                                    hbond_pair[0][2], hbond_pair[0][1], resnum_1, resname_1, # donor_info
                                    hbond_pair[2][1], hbond_pair[2][0], # donor_related_atom_names
                                    residue_1[hbond_pair[2][1]].get_coord(), # D_coord
                                    residue_1[hbond_pair[2][0]].get_coord(), # H_coord
                                    hbond_pair[1][2], hbond_pair[1][1], resnum_2, resname_2, # acceptor_info
                                    hbond_pair[3][0], hbond_pair[3][1], hbond_pair[3][2], # acceptor_related_atom_names 
                                    residue_2[hbond_pair[3][0]].get_coord(), # A_coord
                                    residue_2[hbond_pair[3][1]].get_coord(), # A1_coord
                                    residue_2[hbond_pair[3][2]].get_coord()  # A2_coord
                                )
                            )

                        else: # res 2 is the donor                         
                            if hbond_pair[2][0] == 'HN': 
                                simple_NH_hbond_count[resnum_2] += 1 
                            elif hbond_pair[2][0] == 'H': 
                                simple_CA_hbond_count[resnum_2] += 1
                            else: 
                                simple_SC_hbond_count[resnum_2] += 1

                            if hbond_pair[3][0] == 'O': 
                                simple_CO_hbond_count[resnum_1] += 1
                            else: 
                                simple_SC_hbond_count[resnum_1] += 1

                            hbonds_pairs_in_model.append(
                                hbond_tuple(
                                    hbond_pair[0][2], hbond_pair[0][1], resnum_2, resname_2, # donor_info
                                    hbond_pair[2][1], hbond_pair[2][0], # donor_related_atom_names
                                    residue_2[hbond_pair[2][1]].get_coord(), # D_coord
                                    residue_2[hbond_pair[2][0]].get_coord(), # H_coord
                                    hbond_pair[1][2], hbond_pair[1][1], resnum_1, resname_1, # acceptor_info
                                    hbond_pair[3][0], hbond_pair[3][1], hbond_pair[3][2], # acceptor_related_atom_names
                                    residue_1[hbond_pair[3][0]].get_coord(), # A_coord
                                    residue_1[hbond_pair[3][1]].get_coord(), # A1_coord
                                    residue_1[hbond_pair[3][2]].get_coord()  # A2_coord
                                )
                            )

                    except KeyError:
                        continue

    return hbonds_pairs_in_model, np.array([simple_NH_hbond_count, simple_CO_hbond_count, 
                                            simple_SC_hbond_count, simple_CA_hbond_count]) 


#### <----- ----- ----- ----- burial parameters ----- ----- ----- -----> ####
'''
CB raw burial
Raw count of buried heavy atoms in the CB hemisphere defined by Aashish
Use HA2 burial for GLY's CB burial, hemisphere radius: 5 A 

NH raw burial 
0 for PRO and for C-ter residues 

All parameters here are taken from Aashish's PNAS paper. 
''' 
cb_radius    = 11
cb_ignoreN   = 5
cb_coneangle = 90

nh_radius    = 5
nh_ignoreN   = 0
nh_coneangle = 90

co_radius    = 5
co_ignoreN   = 0
co_coneangle = 90


#### <----- ----- ----- ----- ----- ----- ----- -----> ####
output_head_string  = '#Index chain_id chain_res_id restype secseq ' 
output_head_string += 'CB_burial NH_burial CO_burial SASA ' 
output_head_string += 'NH_hbond CO_hbond SC_hbond CA_hbond '
output_head_string += 'CAx CAy CAz CBx CBy CBz Nx Ny Nz Ox Oy Oz SCMx SCMy SCMz' 


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Prepare input file',
        usage ='use "%(prog)s --help" for more information')

    parser.add_argument('pdb',    help='[required] input .pdb file')
    parser.add_argument('--model',default=0, help='Choose only a specific model in the .pdb')

    parser.add_argument('--membrane-thickness', default=30.0, type=float,
	help='Thickness of membrane bilayer, if the protein is a membrane protein.')

    args = parser.parse_args()

    structure = Bio.PDB.PDBParser(QUIET=True).get_structure('protein',args.pdb)
    model     = structure[args.model] if args.model != None else structure[0]
    nres      = len([res for res in model.get_residues() if res.get_id()[0] == ' '])

    ''' 
    = = = = = = = = = = = = = = = = = = = = 
    coordinates 
    = = = = = = = = = = = = = = = = = = = = 
    '''
    CA_coords     = np.zeros((nres,3))
    CB_coords     = np.zeros((nres,3))
    N_coords      = np.zeros((nres,3))
    O_coords      = np.zeros((nres,3))
    SC_CoM_coords = np.zeros((nres,3))  

    restypes       = ['' for i in range(nres)]
    chain_ids      = ['' for i in range(nres)]
    chain_res_ids  = np.zeros(nres,dtype=int) 

    for resnum, residue in enumerate(model.get_residues()):
        CA_coords[resnum] = residue['CA'].get_coord()
        N_coords[resnum]  = residue['N'].get_coord()
        O_coords[resnum]  = residue['O'].get_coord()
        if residue.get_resname() != 'GLY':
            CB_coords[resnum] = residue['CB'].get_coord() 
        else:
            try:
                CB_coords[resnum] = residue['HA2'].get_coord()
            except KeyError:
                CB_coords[resnum] = residue['CA'].get_coord()
        SC_CoM_coords[resnum] = get_sidechain_com(residue) 
            
        restypes[resnum]  = residue.get_resname() # 3-letter name 
        chain_ids[resnum] = residue.get_full_id()[2]  # get_full_id ("1abc", 0, "A", (" ", 10, "A")) structure, model, chain, residue 
        chain_res_ids[resnum] = residue.get_full_id()[3][1]

    ''' 
    = = = = = = = = = = = = = = = = = = = = 
    CB, NH, CO heavy-atom burial 
    = = = = = = = = = = = = = = = = = = = = 
    '''
    # heavyatom_list includes non-CB heavy atoms from sidechain, 
    # different from Aashish's definition which only has [CB, CA, C, N, O]
    heavyatom_list = []
    for resnum, residue in enumerate(model.get_residues()):
        for atom in residue.get_atom():
            if atom.element != 'H': 
                heavyatom_list.append(atom) 

    CB_rawburial = np.zeros(nres,dtype=int) 
    NH_rawburial = np.zeros(nres,dtype=int)
    CO_rawburial = np.zeros(nres,dtype=int)

    for resnum, residue in enumerate(model.get_residues()):
        resname = residue.get_resname()
        res_chain_id = residue.get_full_id()[2]
        chain_res_id = residue.get_full_id()[3][1]
        
        cb_atom   = residue['CB'] if resname != 'GLY' else residue['HA2']
        cb_radius = 11 if resname != 'GLY' else 5 

        for atom in heavyatom_list:
            atom_chain_id     = atom.get_parent().get_full_id()[2]
            atom_chain_res_id = atom.get_parent().get_full_id()[3][1]
            
            # CB burial 
            if (( res_chain_id!=atom_chain_id) or  # neighbor check
                ((res_chain_id==atom_chain_id) and (math.fabs(chain_res_id-atom_chain_res_id)>cb_ignoreN))): 
                if cb_atom - atom < cb_radius:    # distance check
                    if np.rad2deg(
                        calc_angle(residue['CA'].get_vector(),
                                   cb_atom.get_vector(),
                                   atom.get_vector())) > cb_coneangle: # cone-angle check
                        CB_rawburial[resnum] += 1
            
            # CO burial 
            if (( res_chain_id!=atom_chain_id) or 
                ((res_chain_id==atom_chain_id) and (math.fabs(chain_res_id-atom_chain_res_id)>co_ignoreN))):
                try:
                    o_atom = residue['O']
                    if o_atom - atom < co_radius: 
                        if np.rad2deg(
                            calc_angle(residue['C'].get_vector(),
                                   o_atom.get_vector(),
                                   atom.get_vector())) > co_coneangle:
                            CO_rawburial[resnum] += 1
                except KeyError:
                    CO_rawburial[resnum] = 0
            
            # NH burial 
            try:
                nh_atom = residue['HN']
                if (( res_chain_id!=atom_chain_id) or 
                    ((res_chain_id==atom_chain_id) and (math.fabs(chain_res_id-atom_chain_res_id)>nh_ignoreN))):
                    if nh_atom - atom < nh_radius: 
                        if np.rad2deg(
                            calc_angle(residue['N'].get_vector(),
                                       nh_atom.get_vector(),
                                       atom.get_vector())) > nh_coneangle:
                            NH_rawburial[resnum] += 1
            except KeyError:
                continue 

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
    #dssp = Bio.PDB.DSSP(model,args.pdb)
    #secseq = ['N' for i in range(nres)]
    #sasa   = np.zeros(nres,dtype=int)
    #relative_sasa = np.zeros(nres,dtype='f4') 
    import mdtraj as md 
    secseq = md.compute_dssp(md.load(args.pdb)[args.model], simplified=False)[0]
    secseq = np.where(secseq == ' ', 'N', secseq)
    sasa = md.shrake_rupley(md.load(args.pdb)[args.model], n_sphere_points=960, mode='residue')[0]

    #for i, res in enumerate(list(dssp)):
    #    if res[1] != '-':
    #        secseq[i] = res[1]
    #    sasa[i] = res[2]
    #    relative_sasa[i] = round(res[3],3)

    '''
    = = = = = = = = = = = = = = = = = = = = 
    Hbond 
    = = = = = = = = = = = = = = = = = = = = 
    '''
    hbonds_pairs_in_model, simple_counts = find_hbond_in_model(model) 
    # [simple_NH_hbond_count, simple_CO_hbond_count, simple_SC_hbond_count, simple_CA_hbond_count]

    #frozen = jp.encode(hbonds_pairs_in_model)
    #with open(args.pdb+'.hbond_pairs.pkl','w') as f:
    #    json.dump(frozen,f)

    lb = tb.open_file(bsnm(args.pdb)+'.hbond_pairs.h5',mode='w')
    lb.create_array(lb.root, 'donor_chain_ids',           np.array([p.donor_chain_id           for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'donor_res_index_in_chains', np.array([p.donor_res_index_in_chain for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'donor_res_index_in_models', np.array([p.donor_res_index_in_model for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'donor_resnames',            np.array([p.donor_resname            for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'D_names',  np.array([p.D_name  for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'H_names',  np.array([p.H_name  for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'D_coords', np.array([p.D_coord for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'H_coords', np.array([p.H_coord for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'acceptor_chain_ids',           np.array([p.acceptor_chain_id           for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'acceptor_res_index_in_chains', np.array([p.acceptor_res_index_in_chain for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'acceptor_res_index_in_models', np.array([p.acceptor_res_index_in_model for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'acceptor_resnames',            np.array([p.acceptor_resname            for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A_names',   np.array([p.A_name   for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A1_names',  np.array([p.A1_name  for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A2_names',  np.array([p.A2_name  for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A_coords',  np.array([p.A_coord  for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A1_coords', np.array([p.A1_coord for p in hbonds_pairs_in_model]))
    lb.create_array(lb.root, 'A2_coords', np.array([p.A2_coord for p in hbonds_pairs_in_model]))
    lb.close()

    '''
    = = = = = = = = = = = = = = = = = = = = 
    output
    = = = = = = = = = = = = = = = = = = = = 
    '''
    with open(bsnm(args.pdb) + '.residue_properties.dat','w') as f:
        print >> f, '#thicnkess %.1f' % args.membrane_thickness
        print >> f, output_head_string
        for resnum in range(nres):
            print >> f, '%4i %s %4i %3s %s %3i %3i %3i %6.3f %i %i %i %i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f' % \
                (resnum, chain_ids[resnum], chain_res_ids[resnum], restypes[resnum], secseq[resnum], 
                CB_rawburial[resnum], NH_rawburial[resnum], CO_rawburial[resnum], sasa[resnum],
                simple_counts[0][resnum], simple_counts[1][resnum], simple_counts[2][resnum], simple_counts[3][resnum],
                CA_coords[resnum][0], CA_coords[resnum][1], CA_coords[resnum][2],
                CB_coords[resnum][0], CB_coords[resnum][1], CB_coords[resnum][2],
                N_coords[resnum][0], N_coords[resnum][1], N_coords[resnum][2],
                O_coords[resnum][0], O_coords[resnum][1], O_coords[resnum][2],
                SC_CoM_coords[resnum][0], SC_CoM_coords[resnum][1], SC_CoM_coords[resnum][2])  


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    end_time   = time.time()
    print
    print 'Total running time:'
    print "--- %s seconds ---" % round((end_time - start_time),3)
