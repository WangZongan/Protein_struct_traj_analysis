##############################################################################
# modeller_make_alpha-helix.py                                              
# generate alpha-helix from a given fasta sequence using modeller 
#
# use: python modeller_make_alpha-helix.py name.pdb fasta_seq start end 
#############################################################################

__author__  = 'Wang Zongan'
__version__ = '2015-03-10'

import sys
from modeller import *
from modeller.optimizers import conjugate_gradients
def make_alpha(pdbname,fasta_seq,start,end):
    # Set up environment
    e = environ()
    e.libs.topology.read('/home/zonganw/modeller/modeller-9.14/modlib/top_heav.lib')
    e.libs.parameters.read('/home/zonganw/modeller/modeller-9.14/modlib/par.lib')

    # Build an extended chain model from primary sequence, and write it out
    m = model(e)
    m.build_sequence(fasta_seq)
    #m.write(file='extended-chain.pdb')

    # Make stereochemical restraints on all atoms
    allatoms = selection(m)
    m.restraints.make(allatoms, restraint_type='STEREO', spline_on_site=False)
    m.restraints.add (secondary_structure.alpha(m.residue_range(start,end))  )

    # Get an optimized structure with CG, and write it out
    cg = conjugate_gradients()
    cg.optimize(allatoms, max_iterations=100000)
    m.orient()

    # add pdb remarks 
    m.remark = """REMARK    built from modeller""" 
    outputname = pdbname[:-4] + '_' + str(start) + '_' + str(end) + '.pdb'
    m.write(file=outputname)

def main():
    pdbname,fasta_seq,start,end = sys.argv[1:]
    make_alpha(pdbname,fasta_seq,start,end)

if __name__ == '__main__':
    main()
