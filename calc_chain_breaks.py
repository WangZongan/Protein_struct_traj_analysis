# coding: utf-8
'''
Calculate chain breaks in model.
'''
import re,gc 
import sys,os,time,math,string 
import cPickle as cp
from glob import glob
import collections
import pandas as pd
import numpy as np
np.set_printoptions(precision=3,suppress=True)

__author__  = 'Wang Zongan'
__version__ = '2016-04-21'


def calculate_chain_break(pdb):
    import prody
    structure = prody.parsePDB(pdb)
    sel_str = 'name N CA C'
    coords = structure.select(sel_str).getCoords()
    print len(coords)/3., 'residues found.'
    bond_lengths = np.sqrt(np.sum(np.diff(coords,axis=0)**2,axis=-1))

    if bond_lengths.max() > 2.:
        breaks = bond_lengths>2.
        print "WARNING: %i separate chains found" % (1+breaks.sum())
        with open(pdb+'.chain_break.chain_first_residues.dat','w') as f:
            outputstring = ''
            for br in list(breaks.nonzero()[0]/3):
                outputstring += '%i,' % (br+1)
            outputstring = outputstring.strip(',')
            f.write(outputstring)

        with open(pdb+'.chain_break.chain_first_residues_config.dat','w') as f:
            outputstring = ''
            for br in list(breaks.nonzero()[0]/3):
                outputstring += '--chain-first-residue=%i ' % (br+1)
            f.write(outputstring)
        
        with open(pdb+'.chain_break.hbond_exclude_residues.dat','w') as f:
            outputstring=','.join(str(list(breaks.nonzero()[0]/3))[1:-1].split(', '))
            for br in list(breaks.nonzero()[0]/3):
                outputstring += ',%i' % (br-1)
                outputstring += ',%i' % (br+1)
            f.write(outputstring)
        
        with open(pdb+'.chain_break.chain_restraint_groups_config.dat','w') as f:
            br_list = [0]
            for br in list(breaks.nonzero()[0]/3):
                br_list.append(br+1)  # all first residues in chain
            br_list.append(int(len(coords)/3.))
            outputstring = ''
            print "Last residue of each chain: %s " % br_list[1:] 
            for ibr in range(len(br_list)-1):
                outputstring += '--restraint-group=%i-%i ' % (br_list[ibr],br_list[ibr+1]-1)
            f.write(outputstring)

def main():
    pdb = sys.argv[1]
    calculate_chain_break(pdb)
    
if __name__ == "__main__":
    main()

