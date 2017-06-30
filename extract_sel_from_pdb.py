#!/usr/bin/env python
__author__  = 'Wang Zongan'
__version__ = '2016.09.20'

import os
import re
import sys
import string
import numpy as np
from Bio import PDB


class ChainExtracter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser(QUIET=True)
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, out_path, chain_letters, overwrite=False):
        """ 
        Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        chain_letters = [chain.upper() for chain in chain_letters]
        pdb_fn = os.path.split(pdb_path)[1]
        
        print "OUT PATH:",out_path

        # Skip PDB generation if the file already exists
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_fn, out_path))
            return out_path
        print("Extracting chain%s %s from %s..." % (plural, ", ".join(chain_letters), pdb_fn))

        # Get structure, write new file with only given chains
        struct = self.parser.get_structure('protein', pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


class ResidueExtracter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser(QUIET=True)
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, out_path, residue_indices, overwrite=False):
        """ 
        Create a new PDB file containing only the specified residues.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param residue_indices: np array of selected residue indices 
        :param overwrite: write over the output file if it exists
        """
        residue_indices = np.array(residue_indices)
        pdb_fn = os.path.split(pdb_path)[1]
        
        print "OUT PATH:",out_path

        # Skip PDB generation if the file already exists
        plural = "s" if (len(residue_indices) > 1) else ""  # for printing
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Residue%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(residue_indices), pdb_fn, out_path))
            return out_path
        print("Extracting %i residue%s \n%s from %s..." % (len(residue_indices), plural, ", ".join(residue_indices.astype(str)), pdb_fn))

        # Get structure, write new file with only given chains
        struct = self.parser.get_structure('protein', pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectResidues(residue_indices))

        return out_path


class SelectResidues(PDB.Select):
    """ Only accept the specified residues when saving. """
    def __init__(self, residue_indices):
        self.residue_indices = residue_indices

    def accept_residue(self, residue):
        return (residue.get_id()[1]-1 in self.residue_indices)


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


def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage ='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter,
        description = 'Extract selected chain(s) or residue(s) from given pdb file.')
    parser.add_argument('pdbpath', help = textwrap.dedent('''[required] Path to pdb trajectories.'''))
    parser.add_argument('--residues', default=None, type=parse_segments,
        help = textwrap.dedent('''
            If provided, only the selected residues will be extracted to a new pdb file; otherwise, all residues are selected.
            The selection if of the form 10-30,50-60. Caution: the pdb file better not have HETATM or water molecures.'''))
    parser.add_argument('--chains', default=[], action='append', type=str, 
        help = textwrap.dedent('''
            If provided, only the selected chains will be extracted to a new pdb file; otherwise all chains are selected. ''')) 
    parser.add_argument('--output-path', default=None, type=str, help = textwrap.dedent('''Output path.''')) 
    args = parser.parse_args()


    if args.chains:
        extracter = ChainExtracter(args.output_path)
        extracter.make_pdb(args.pdbpath, args.output_path, args.chains, overwrite=True)

    if args.residues is not None:
        extracter = ResidueExtracter(args.output_path)
        extracter.make_pdb(args.pdbpath, args.output_path, args.residues, overwrite=True)


if __name__ == "__main__":
    main()
