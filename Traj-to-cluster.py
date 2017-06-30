#!/usr/bin/env python
'''
Given pdb trajectories, perform clustering analysis.
'''
__author__  = 'Wang Zongan'
__version__ = '2016-12-11'

import os
import sys
import time
import string
import numpy as np
import pandas as pd
import cPickle as cp
import mdtraj as md
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def featurize_trajectories(coords, featurizer):
    '''
    Input
    coords : list of 'MDTrajDataset' object

    Output 
    features : list of arrays, length n_trajs, each of shape (n_samples, n_features)
    '''
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
	
	
def decompose_features(features, decomposer, n_components=None, lag_time=1):
    '''
    Decomposing features is a way to reduce the dimension of the features. 

    Each of the components is a eigenvector of the feature space, dimension: (n_features,) 

    The old features are transformed to the new feature space. 
    
    Consider one sample, which is vectorized to (n_features,).T, 
    apply the transform matrix, which is in the shape (n_components, n_features), 
    we will get its projection onto the new space (n_components,). 

    --------------------------------------------------------------------------------------------------------------------------------------
    Input
    features         : array-like, length n_trajs, each of shape (n_samples, n_features)
	
    Output
    features_new     : array-like, length n_trajs, each of shape (n_samples, n_components) ((n_samples, n_samples) if n_components = None) 

    dcmp.components_ : shape (n_components, n_features), ((n_samples, n_features) if n_components = None)
        PCA  : Principal axes in feature space, representing the directions of maximum variance in the data.
        tICA : Components with maximum autocorrelation. 
    '''
    if decomposer == 'PCA':
	from msmbuilder.decomposition import PCA
	dcmp = PCA(n_components=n_components)
    elif decomposer == 'tICA':
	from msmbuilder.decomposition import tICA
	dcmp = tICA(n_components=n_components, lag_time=lag_time)
    features_new = dcmp.fit_transform(features)
    return features_new, dcmp.components_

	
def cluster_features(features, clusterer, n_clusters=8):
    '''
    Input
    features : list of arrays, length n_trajs, each of shape (n_samples, n_features)
	
    Output
    clst : msmbuilder.cluster object, with attributes
        cluster_centers_ : (n_clusters, n_features)
	labels_	         : list of arrays, each of shape (n_samples, )
    '''
    if clusterer == 'KMeans':
	from msmbuilder.cluster import KMeans
	clst = KMeans(n_clusters=n_clusters)
    elif clusterer == 'KCenters':
	from msmbuilder.cluster import KCenters
	clst = KCenters(n_clusters=n_clusters)
    elif clusterer == 'KMedoids':
	from msmbuilder.cluster import KMedoids
	clst = KMedoids(n_clusters=n_clusters)
    elif clusterer == 'MiniBatchKMeans':
	from msmbuilder.cluster import MiniBatchKMeans
	clst = MiniBatchKMeans(n_clusters=n_clusters)
    elif clusterer == 'MiniBatchKMedoids':
	from msmbuilder.cluster import MiniBatchKMedoids
        clst = MiniBatchKMedoids(n_clusters=n_clusters)	
    clusters = clst.fit_transform(features)
    return clst


def cluster_centroid_size(labels):
    from collections import Counter
    return Counter(np.array(labels)).most_common()

	
def find_cluster_center_indices(features, clst):
    '''
    Input
    features : list of arrays, length n_trajs, each of shape (n_samples, n_features)
    clst     : msmbuilder.cluster object, with attributes
        cluster_centers_ : (n_clusters, n_features)
	labels_		 : list of arrays, each of shape (n_samples, )
	
    Output
    cci : indices of structures among samples that are closest to the cluster centers, if they are not in the samples.
    '''
    from scipy.spatial.distance import cdist
    cat_features = np.concatenate(features)     # (n_trajectories * n_samples, n_features)
    true_cc      = clst.cluster_centers_        # (n_clusters, n_features)
    cat_features = cdist(cat_features, true_cc) # (n_trajectories * n_samples, n_clusters)
    cci = [np.where(cat_features[:,i] == cat_features[:,i].min())[0][0] for i in range(len(true_cc))]
    return np.array(cci)


def plot_cluster(X, Y, cluster_center_indices, figname, x_label, y_label, xmin, xmax, ymin, ymax, c_map='winter', cc_color='red'):
    fig = figure(figsize=(10,10))
    hexbin(X, Y, bins='log', mincnt=1, cmap=c_map) 
    # mincnt: [ None | a positive integer ] if not None, only display cells with more than mincnt number of points in the cell
    scatter(X[cluster_center_indices], Y[cluster_center_indices], s=100, c=cc_color)
    xlabel(x_label, fontsize=20)
    ylabel(y_label, fontsize=20)
    xlim(xmin, xmax) 
    ylim(ymin, ymax) 
    tick_params(axis='both', which='major', labelsize=20)
    grid()
    tight_layout()
    savefig(figname, format='png')

	
def ceil(number, decimal=1):
    return np.ceil(number*np.power(0.1,decimal))*np.power(10,decimal)
	
def floor(number, decimal=1):
    return np.floor(number*np.power(0.1,decimal))*np.power(10,decimal) 


def toNumpy32(dataset):
    '''
    dataset : array-like
    '''
    return [ds.astype(np.float32) for ds in dataset] 


bsnm = lambda fpath: os.path.splitext(os.path.basename(fpath))[0]
	

def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage = textwrap.dedent('''Use "python %(prog)s -h" for more information.'''),
        formatter_class=argparse.RawTextHelpFormatter, 
        description = textwrap.dedent('''\
            First, this program employs msmbuilder to featurize given pdb trajectories into vectorizable space.
	    Second, the vector space is decompose by tICA or PCA to further reduce the dimension. 
            Third, clustering is performed so that each structure in the trajectories is labeled by an index. 
        
            Example:
            $ python Traj-to-cluster.py     \n
                path_to_pdb_trajectories/   \n
                --featurizer=DRIDFeaturizer \n
                --decomposer=PCA            \n
                --decomposer-n-components=5 \n
                --clusterer=KCenters        \n
                --n-clusters=5 '''))
    parser.add_argument('pdbpath', help = textwrap.dedent('''[required] Path to pdb trajectories.'''))
    parser.add_argument('--lag-time', default=1, type=int, help = textwrap.dedent('''Lag time of the model. Default value = 1.'''))
    parser.add_argument('--featurizer', default=None, type=str, 
        help = textwrap.dedent('''\
            Featurizer at your choice. Available featurizers are (select them by name): 
            (1) RMSDFeaturizer;
            (2) DihedralFeaturizer, only phi and psi angles;
            (3) DRIDFeaturizer (DRID, Distribution of Reciprocal of Interatomic Distances);
            (4) ContactFeaturizer, CA contact. 	
            Note: user must choose a featurization method. Choose by name. '''))
    parser.add_argument('--decomposer', default=None, type=str, 
	help = textwrap.dedent('''\
            Decomposer at your choice. Available decomposers are: 
            (1) tICA;
            (2) PCA. 
            Note: selection of decomposer is not necessary but recommended.
            If not provided, program will ignore this step and cluster directly on raw features. '''))
    parser.add_argument('--decomposer-n-components', default=None, type=int,  
	help = textwrap.dedent('''Number of components to keep. if n_components is not set all components are kept.'''))
    parser.add_argument('--clusterer', default=None, type=str,
	help = textwrap.dedent('''\
            Clustering method at your choice. Available clusterer are: 
            (1) KMeans;
            (2) KCenters;
            (3) KMedoids;
            (4) MiniBatchKMeans;
            (5) MiniBatchKMedoids.
            Note: user must choose a clusering method. '''))
    parser.add_argument('--n-clusters', default=5, type=int, 
	help = textwrap.dedent('''The number of clusters to form as well as the number of centroids to generate.'''))
    parser.add_argument('--reference-model', default=[], action='append', type=str, 
        help = textwrap.dedent(''' Reference models used to calculate RMSD. '''))
    parser.add_argument('--output', type=str, default='output', help = textwrap.dedent('''Output file name.'''))
    args = parser.parse_args()

    from msmbuilder.dataset import dataset
    coords = dataset(os.path.join(args.pdbpath, '*.pdb')) # coords: 'MDTrajDataset' object 
    print '%i trajectories found. \n' % len(coords)

    ## featurize
    features = featurize_trajectories(coords, args.featurizer)
    print "%s selected" % args.featurizer
    print "features: (n_samples, n_features) = (%i, %i) for each trajectory \n" % (features[0].shape[0], features[0].shape[1])
    with open(args.output+'.features.%s.pkl' % args.featurizer, 'w') as f:
        cp.dump(toNumpy32(features), f, -1)
    sys.stdout.flush()
	
    ## decompose
    if args.decomposer == None:
	print "No decomposer is selected! Program will directly cluster the raw coordinates. \n"
    else:
        features, components = decompose_features(features, args.decomposer, n_components=args.decomposer_n_components, lag_time=args.lag_time)
	print "%s selected" % args.decomposer
	print "features: (n_samples, n_features) = (%i, %i) for each trajectory \n" % (features[0].shape[0], features[0].shape[1]) 
		
    ## clustering
    clst = cluster_features(features, args.clusterer, n_clusters=args.n_clusters)
    cci  = find_cluster_center_indices(features, clst)
    print "%s selected" % args.clusterer
    print "Cluster center indices: %s \n" % cci 

    cat_features = np.concatenate(features)
    cat_labels   = np.concatenate(clst.labels_)

    ## reference pdb
    if args.reference_model != None:
        import mdtraj as md
        from msmbuilder.featurizer import RMSDFeaturizer
        rmsd_to_ref = [] 
        for ref in args.reference_model:
            print "\nCompute RMSD to the reference models : %s." % ref
            ref_traj = md.load(ref)
            print "N atoms %i" % ref_traj.n_atoms
            rmsd_to_ref.append(np.concatenate(RMSDFeaturizer(ref_traj[0]).fit_transform(coords))[:,0])
        rmsd_to_ref = np.array(rmsd_to_ref)*10 # unit : A

        for i in range(len(args.reference_model)):
            with open(args.output+'.ref_RMSD_%s.dat'%bsnm(args.reference_model[i]), 'w') as f:
                for e in rmsd_to_ref[i]:
                    print >> f, '%.3f' % e
        
        plot_cluster(X=rmsd_to_ref[0], Y=rmsd_to_ref[1], cluster_center_indices=cci,
                figname=args.output+'.ref_RMSD.png',
                x_label='state 1', y_label='state 2',
                xmin=0, xmax=ceil(rmsd_to_ref.max()),
                ymin=0, ymax=ceil(rmsd_to_ref.max()),
                c_map='winter', cc_color='red')

    ## output 
    with open(args.output+'.labels.pkl', 'w') as f:
        cp.dump(cat_labels, f, -1)

    with open(args.output+'.components.pkl', 'w') as f:
        cp.dump(toNumpy32(components), f, -1)

    with open(args.output+'.cluster_center_idx.dat', 'w') as f:
        for i in range(args.n_clusters):
            print >> f, '%6i' % cci[i]
	
    if args.decomposer == 'tICA':
	plot_cluster(X=cat_features[:,0], Y=cat_features[:,1], cluster_center_indices=cci, 
	    figname=args.output+'.tICA_1st_2nd.png', 
	    x_label='tIC 1', y_label='tIC 2', 
	    xmin=floor(cat_features[:,0].min()), xmax=ceil(cat_features[:,0].max()),
	    ymin=floor(cat_features[:,1].min()), ymax=ceil(cat_features[:,1].max()),
	    c_map='winter', cc_color='red')
    elif args.decomposer == 'PCA':
	plot_cluster(X=cat_features[:,0], Y=cat_features[:,1], cluster_center_indices=cci, 
	    figname=args.output+'.PCA_1st_2nd.png', 
	    x_label='PC 1', y_label='PC 2',
	    xmin=floor(cat_features[:,0].min()), xmax=ceil(cat_features[:,0].max()),
	    ymin=floor(cat_features[:,1].min()), ymax=ceil(cat_features[:,1].max()),
	    c_map='winter', cc_color='red')


def sec_to_hr_min_sec(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


if __name__ == '__main__':
    sta = time.time()
    main()
    print '\nrunning time: %s\n' % sec_to_hr_min_sec(time.time() - sta)


