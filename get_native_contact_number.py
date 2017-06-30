#!/usr/bin/env python
'''
'''
__author__  = 'Wang Zongan'
__version__ = '2016-10-04'

import os
import sys
import string
import numpy as np
import pandas as pd
import cPickle as cp
import mdtraj as md
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def featurize_trajectories(coords, featurizer):
    if featurizer == 'RMSDFeaturizer':
	from msmbuilder.featurizer import RMSDFeaturizer
	feat = RMSDFeaturizer(reference_traj=coords[0])		
    elif featurizer == 'DRIDFeaturizer':
        from msmbuilder.featurizer import DRIDFeaturizer
        feat = DRIDFeaturizer()	
    elif featurizer == 'ContactFeaturizer':
        from msmbuilder.featurizer import ContactFeaturizer
	feat = ContactFeaturizer(scheme='ca')	
    elif featurizer == 'DihedralFeaturizer':
	from msmbuilder.featurizer import DihedralFeaturizer
	feat = DihedralFeaturizer(types=['phi', 'psi'])
    return feat.fit_transform(coords)	
	

def get_basename_no_ext(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage = textwrap.dedent('''Use "python %(prog)s -h" for more information.'''),
        formatter_class=argparse.RawTextHelpFormatter) 
    parser.add_argument('pdbpath', help = textwrap.dedent('''[required] Path to pdb trajectories.'''))
    parser.add_argument('target', help = textwrap.dedent('''[required] Path to target pdb.
        Note: The target pdb should have the same number of atoms in structure with that in pdb trajectories. '''))
    args = parser.parse_args()

    from msmbuilder.dataset import dataset
    coords = dataset(args.pdbpath) 
    print '%i trajectories found. ' % len(coords)
	
    ## featurize
    features = featurize_trajectories(coords, 'ContactFeaturizer')
    #print "features: (n_samples, n_features) = (%i, %i) for each trajectory \n" % (features[0].shape[0], features[0].shape[1])
	
    import mdtraj as md
    target = md.load(args.target)
		
    native_contact_dists, native_contact_pairs = md.compute_contacts(target, scheme='ca')
    native_contact_pairs = native_contact_pairs[np.where(native_contact_dists[0]<=0.75)]
    n_native_contact = len(native_contact_pairs)
    print "Target structure has %i pairs of CA-CA contact in total. \n" % n_native_contact
		
    from msmbuilder.featurizer import ContactFeaturizer
    native_contact_to_target = np.concatenate(ContactFeaturizer(contacts=native_contact_pairs,scheme='ca').fit_transform(coords)) # (n_samples, n_pairs)
    native_contact_to_target = np.select([native_contact_to_target<=0.75,native_contact_to_target>0.75],[1,0]) 
    native_contact_to_target = np.sum(native_contact_to_target,axis=1)
    
    with open('%s.%s.number_native_contact.dat' % (get_basename_no_ext(args.target), get_basename_no_ext(args.pdbpath)), 'w') as f:
        for e in native_contact_to_target:
            print >> f, '%i %i %.3f' % (n_native_contact, e, e*1./n_native_contact)

if __name__ == '__main__':
    main()
