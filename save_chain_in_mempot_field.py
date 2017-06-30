#!/usr/bin/env python
'''
'''

__author__ = 'Wang Zongan'
__version__ = '2017.04.25'

import sys
import Bio.PDB
from Bio.PDB.PDBIO import Select

import time
import StringIO
import os,math
import numpy as np


class ChainInMemPotField(Select):
    def __init__(self, thickness, field):
        self.thickness = thickness
        self.field     = field

    def accept_chain(self, chain):
        CAz = np.array([res['CA'].get_coord()[2] for res in chain])
        CAz_included = len(np.where(
                            (CAz <=  self.thickness/2. + self.field) & 
                            (CAz >= -self.thickness/2. - self.field))[0])
        if CAz_included <= 5: 
            return 0
        else:
            print 'Save chain %s' % chain.get_id()
            return 1


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Program for only saving chains in the membrane potential field.',
        usage ='use "python %(prog)s --help or -h" for more information')
    parser.add_argument('inputpdb', help='[required] input pdb file')
    parser.add_argument('--thickness', type=float, help='hydrophobic thickness of lipid bilayer')
    parser.add_argument('--field', type=float, help='potential field')
    args = parser.parse_args()

    protein = Bio.PDB.PDBParser(QUIET=True).get_structure('x',args.inputpdb)

    io = Bio.PDB.PDBIO()
    io.set_structure(protein)
    io.save(bsnm(args.inputpdb) + '_InMemPotField.pdb', select=ChainInMemPotField(args.thickness, args.field)) 


if __name__ == '__main__':
    main()

