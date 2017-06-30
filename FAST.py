#!/usr/bin/env python
'''
Given pdb trajectories, perform FAST analysis.
'''
__author__  = 'Wang Zongan'
__version__ = '2016-09-05'

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
	
	
def decompose_features(features, decomposer, n_components=None, lag_time=1):
    '''
    Input
    features : list of arrays, length n_trajs, each of shape (n_samples, n_features)
	
    Output
    features_new : list of arrays, length n_trajs, each of shape (n_samples, n_features_new) 
    '''
    if decomposer == 'PCA':
	from msmbuilder.decomposition import PCA
	dcmp = PCA(n_components=n_components)
    elif decomposer == 'tICA':
	from msmbuilder.decomposition import tICA
	dcmp = tICA(n_components=n_components, lag_time=lag_time)
    return dcmp.fit_transform(features)

	
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
    cci : indices of structures among samples that are closest to the cluster centers, if they are not in the samples
    '''
    from scipy.spatial.distance import cdist
    cat_features = np.concatenate(features)     # (n_trajectories * n_samples, n_features)
    true_cc      = clst.cluster_centers_        # (n_clusters, n_features)
    cat_features = cdist(cat_features, true_cc) # (n_trajectories * n_samples, n_clusters)
    cci = [np.where(cat_features[:,i] == cat_features[:,i].min())[0][0] for i in range(len(true_cc))]
    return np.array(cci)


def find_transitions(labels, i, j, lag_time=1):
    '''
    Find the number of transitions from state i to state j in labels.
    '''
    assert i in labels
    assert j in labels
    i_idx = np.where(labels == i)[0]
    i_idx_n = i_idx + lag_time
    i_idx_n = i_idx_n[np.where(i_idx_n < len(labels))]
    i_to_j_idx = i_idx[np.where(labels[i_idx_n] == j)]
    return len(i_to_j_idx)


def calc_transition_count_mat(labels, n_clusters):
    print labels
    transition_count_mat = np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        transition_count_mat[i][i] = len(np.where(labels==i)[0])
        for j in range(i+1, n_clusters):
            transition_count_mat[i][j] = find_transitions(labels,i,j,1)
            transition_count_mat[j][i] = find_transitions(labels,j,i,1)
    return transition_count_mat


def build_msm(labels, lag_time=1, n_timescales=None, prior_counts=0):
    '''
    Input
    labels : list of arrays, each of shape (n_samples, )
	
    Output
    msm : msmbuilder.msm.msm object, with attributes
        n_states_    : The number of states in the model, determined by labels 
	mapping_
	countsmat_   : (n_states_, n_states_), Number of transition counts between states. 
	transmat_    : (n_states_, n_states_), Maximum likelihood estimate of the reversible transition matrix.
	populations_ : (n_states_,), The equilibrium population (stationary eigenvector) of transmat_.
	timescales_
    '''
    from msmbuilder.msm import MarkovStateModel
    msm = MarkovStateModel(lag_time=lag_time, n_timescales=n_timescales, prior_counts=prior_counts)
    msm.fit(labels)
    return msm
	

def calc_FAST_reward_score(physical_scores, cluster_center_indices, transition_count_mat, alpha=1, n_simulations=30, minmax='min'): 
    '''
    physical_scores      : (n_samples,)
    '''
    n_clusters = len(cluster_center_indices) 

    ps_max = physical_scores[cluster_center_indices].max()
    ps_min = physical_scores[cluster_center_indices].min()
    ps_max_min = ps_max - ps_min

    c_max = np.diag(transition_count_mat).max()
    c_min = np.diag(transition_count_mat).min()
    c_max_min = c_max - c_min
    c = np.sum(transition_count_mat, axis=0) - np.diag(transition_count_mat) # total transition from state i to any other state

    if minmax == 'min':
        rewards = np.array([
            ( (ps_max - physical_scores[cluster_center_indices[i]]) / ps_max_min + alpha * (c_max - c[i]) / c_max_min )
            for i in range(n_clusters) ])
    elif minmax == 'max':
        rewards = np.array([
            ( (physical_scores[cluster_center_indices[i]] - ps_min) / ps_max_min + alpha * (c_max - c[i]) / c_max_min )
            for i in range(n_clusters) ])
    else:
        print "minmax must be either min or max!"


    sims = np.around(n_simulations * rewards / np.sum(rewards))
    sims[np.where(sims==0)] = 1 # remove zeros
    if sims.sum() < n_simulations:
        while sims.sum() < n_simulations:
            idx = np.random.randint(n_clusters)
            sims[idx] += 1
    elif sims.sum() > n_simulations:
        while sims.sum() > n_simulations:
            idx = np.random.randint(n_clusters)
            if sims[idx] > 1:
                sims[idx] -= 1
    assert sims.sum() == n_simulations
	
    return rewards, sims, c

	
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
	

def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage = textwrap.dedent('''Use "python %(prog)s -h" for more information.'''),
        formatter_class=argparse.RawTextHelpFormatter, 
        description = textwrap.dedent('''\
            First, this program employs msmbuilder to featurize given pdb trajectories into vectorizable space.
	    Second, the vector space is decompose by tICA or PCA to further reduce the dimension. 
            Third, clustering is performed so that each structure in the trajectories is labeled by an index. 
	    Forth, Marcov State Model, albeit may not be well behaved, is built on the labeled trajectories.
	    Last, FAST reward scores are calculated based on the transition-count matrix and user-chosen physical traits. 
        
            Example:
            $ python FAST.py path_to_pdb_trajectories/ --featurizer=DRIDFeaturizer --decomposer=PCA --decomposer-n-components=5 --clusterer=KCenters --n-clusters=5 --msm-prior-counts=0.2 --physical-trait=target-RMSD --target-pdb=/path_to_target_pdb/target.pdb '''))
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
    parser.add_argument('--msm-n-timescales', default=None, type=int, 
	help = textwrap.dedent('''\
	    The number of dynamical timescales to calculate when diagonalizing the transition matrix. 
	    If not specified, it will compute n_states - 1. '''))
    parser.add_argument('--msm-prior-counts', default=0, type=float,
	help = textwrap.dedent('''\
	    Add a number of 'pseudo counts' to each entry in the counts matrix after ergodic trimming. 
	    When prior_counts == 0 (default), the assigned transition probability between two states 
	    with no observed transitions will be zero, whereas when prior_counts > 0, even this unobserved 
	    transitions will be given nonzero probability. '''))
    parser.add_argument('--physical-trait', default=None, type=str,
	help = textwrap.dedent('''\
            Physical trait used in calculation of FAST reward score. Available choices are: 
            (1) target-RMSD, if chosen, user must supply a target structure; 
            (2) target-native-contact, if chosen, user must supply a target structure; 
            (3) target-tmscore, if chosen, user must supply the data file containing the TM-scores in column;
            (4) potential, target free, if chosen, user must supply the data file containing the potentials in column; 
            Note: user must choose a physical trait. '''))
    parser.add_argument('--target-pdb', default=None, type=str,
	help = textwrap.dedent('''\
            The target pdb structure. 
            Note: The target pdb should have the same number of atoms in structure with that in pdb trajectories. '''))
    parser.add_argument('--initial-pdb', default=None, type=str,
	help = textwrap.dedent('''\
            The initial pdb structure. 
            Note: The initial pdb should have the same number of atoms in structure with that in pdb trajectories. '''))
    parser.add_argument('--potential', default=None, type=str, 
        help = textwrap.dedent('''The potential file corresponding to the pdb trajectories. '''))
    parser.add_argument('--tmscore', default=None, type=str, help = textwrap.dedent('''The TM-score file corresponding to the pdb trajectories. '''))
    parser.add_argument('--fast-n-simulations', default=30, type=int,
        help = textwrap.dedent('''Number of parallel simulations in each round of FAST algorithm. Default value: 30. '''))
    parser.add_argument('--fast-alpha', default=1., type=float, help = textwrap.dedent('''Number of clusters. Default value: 1.0.'''))
    parser.add_argument('--output', type=str, default='output', help = textwrap.dedent('''Output file name.'''))
    args = parser.parse_args()

    from msmbuilder.dataset import dataset
    coords = dataset(os.path.join(args.pdbpath, '*.pdb')) 
    print '%i trajectories found. \n' % len(coords)
	
    ## featurize
    features = featurize_trajectories(coords, args.featurizer)
    print "%s selected" % args.featurizer
    print "features: (n_samples, n_features) = (%i, %i) for each trajectory \n" % (features[0].shape[0], features[0].shape[1])
	
    ## decompose
    if args.decomposer == None:
	print "No decomposer is selected! Program will directly cluster the raw features. \n"
    else:
        features = decompose_features(features, args.decomposer, n_components=args.decomposer_n_components, lag_time=args.lag_time)
	print "%s selected" % args.decomposer
	print "features: (n_samples, n_features) = (%i, %i) for each trajectory \n" % (features[0].shape[0], features[0].shape[1]) 
		
    ## clustering
    clst = cluster_features(features, args.clusterer, n_clusters=args.n_clusters)
    cci = find_cluster_center_indices(features, clst)
    print "%s selected" % args.clusterer
    print "Cluster center indices: %s \n" % cci 
	
    ## build msm 
    #msm = build_msm(clst.labels_, lag_time=args.lag_time, n_timescales=args.msm_n_timescales, prior_counts=args.msm_prior_counts)
    #print msm, '\n'
    #print "Transition count matrix: \n %s \n" % msm.countsmat_
    #print "Relative population of each state: %s \n" % msm.populations_

    ## construct transition count matrix 
    transition_count_mat = calc_transition_count_mat(np.concatenate(clst.labels_), args.n_clusters)
    print 'Transition count matrix: \n', transition_count_mat

    #### calculate FAST reward score
    output_df = pd.DataFrame()
    output_df['idx'] = cci
    output_df['#cluster'] = transition_count_mat.diagonal()

    if args.initial_pdb != None:
        import mdtraj as md 
        initial = md.load(args.initial_pdb)

        from msmbuilder.featurizer import RMSDFeaturizer
        rmsd_to_initial = np.concatenate(RMSDFeaturizer(initial).fit_transform(coords))[:,0]
        
        output_df['iniRMSD'] = rmsd_to_initial[cci]

    if args.target_pdb != None:
        import mdtraj as md
	target = md.load(args.target_pdb)
		
	from msmbuilder.featurizer import RMSDFeaturizer
        rmsd_to_target = np.concatenate(RMSDFeaturizer(target).fit_transform(coords))[:,0]
		
	native_contact_dists, native_contact_pairs = md.compute_contacts(target, scheme='ca')
	native_contact_pairs = native_contact_pairs[np.where(native_contact_dists[0]<=0.75)]
	print "Target structure has %i pairs of CA-CA contact in total. \n" % len(native_contact_pairs)
		
	from msmbuilder.featurizer import ContactFeaturizer
	native_contact_to_target = np.concatenate(ContactFeaturizer(contacts=native_contact_pairs,scheme='ca').fit_transform(coords)) # (n_samples, n_pairs)
	native_contact_to_target = np.select([native_contact_to_target<=0.75,native_contact_to_target>0.75],[1,0]) 
	native_contact_to_target = np.sum(native_contact_to_target,axis=1)

        output_df['tarRMSD'] = rmsd_to_target[cci] 
        output_df['#NativeContact'] = native_contact_to_target[cci]
		
    if args.potential != None:
	potential = np.loadtxt(args.potential)
        output_df['potential'] = potential[cci]
    
    if args.tmscore != None:
	tmscore = np.loadtxt(args.tmscore)
        output_df['tmscore'] = tmscore[cci]
    
    # choose physical trait
    print "%s is selected in FAST \n" % args.physical_trait
    if args.physical_trait == 'target-RMSD':
	if args.target_pdb == None:
	    print "User must provide a target structure! \n"
	rewards, sims, c = calc_FAST_reward_score(rmsd_to_target, cci, transition_count_mat, 
                alpha=args.fast_alpha, n_simulations=args.fast_n_simulations, minmax='min')

    elif args.physical_trait == 'target-native-contact':
	if args.target_pdb == None:
	    print "User must provide a target structure! \n"
	rewards, sims, c = calc_FAST_reward_score(native_contact_to_target, cci, transition_count_mat, 
                alpha=args.fast_alpha, n_simulations=args.fast_n_simulations, minmax='max')

    elif args.physical_trait == 'target-tmscore':
        if args.tmscore == None:
            print "User must provide a TM-score file corresponding to the pdb trajectories! \n"
        rewards, sims, c = calc_FAST_reward_score(tmscore, cci, transition_count_mat, 
                alpha=args.fast_alpha, n_simulations=args.fast_n_simulations, minmax='max')

    elif args.physical_trait == 'potential':
	if args.potential == None:
	    print "User must provide a potential file corresponding to the pdb trajectories! \n"
	rewards, sims, c = calc_FAST_reward_score(potential, cci, transition_count_mat, 
                alpha=args.fast_alpha, n_simulations=args.fast_n_simulations, minmax='min')

    output_df['#Transition'] = c
    output_df['reward'] = rewards
    output_df['#sim'] = sims

    ## output 
    with open(args.output+'.CenterIdx_ClusterSize.dat', 'w') as f:
        for i in range(args.n_clusters):
            print >> f, '%6i %6i' % (cci[i], sims[i])
	
    if args.initial_pdb != None:
        with open(args.output+'.iniRMSD.dat','w') as f:
            for ele in rmsd_to_initial:
                print >> f, '%8.3f' % ele

    if args.target_pdb != None:
	with open(args.output+'.tarRMSD.dat','w') as f:
	    for ele in rmsd_to_target:
		print >> f, '%8.3f' % ele

	with open(args.output+'.tarNativeContact.dat','w') as f:
	    for ele in native_contact_to_target:
		print >> f, '%8.3f' % ele

    with open(args.output+'.dat', 'w') as f:
        print >> f, output_df

    ## plot
    if args.target_pdb != None:
	plot_cluster(X=rmsd_to_target, Y=native_contact_to_target, cluster_center_indices=cci, 
	    figname=args.output+'.tarRMSD_tarNativeContact.png', 
	    x_label='RMSD to target / nm', y_label='# native contact', 
	    xmin=0, xmax=ceil(rmsd_to_target.max(),0),
	    ymin=0, ymax=ceil(native_contact_to_target.max()), 
	    c_map='winter', cc_color='red')
        if args.initial_pdb != None:
            plot_cluster(X=rmsd_to_initial, Y=rmsd_to_target, cluster_center_indices=cci,
                figname=args.output+'.tarRMSD_iniRMSD.png',
	        x_label='RMSD to initial / nm', y_label='RMSD to target / nm', 
	        xmin=0, xmax=ceil(rmsd_to_target.max(),0),
	        ymin=0, ymax=ceil(rmsd_to_initial.max(),0), 
	        c_map='winter', cc_color='red')
        if args.tmscore != None:
	    plot_cluster(X=tmscore, Y=native_contact_to_target, cluster_center_indices=cci, 
		figname=args.output+'.tmscore_tarNativeContact.png', 
		x_label='TM-score to target', y_label='# native contact',
		xmin=0, xmax=1,
		ymin=0, ymax=ceil(native_contact_to_target.max()), 
		c_map='winter', cc_color='red')
            if args.potential != None:
	        plot_cluster(X=tmscore, Y=potential, cluster_center_indices=cci, 
		    figname=args.output+'.tmscore_potential.png', 
		    x_label='TM-score to target', y_label='potential',
		    xmin=0, xmax=1,
		    ymin=floor(potential.min()), ymax=ceil(potential.max()), 
		    c_map='winter', cc_color='red')
        if args.potential != None:
	    plot_cluster(X=rmsd_to_target, Y=potential, cluster_center_indices=cci, 
		figname=args.output+'.tarRMSD_potential.png', 
		x_label='RMSD to target / nm', y_label='potential',
		xmin=0, xmax=ceil(rmsd_to_target.max(),0),
		ymin=floor(potential.min()), ymax=ceil(potential.max()), 
		c_map='winter', cc_color='red')
	
    if args.decomposer == 'tICA':
	cat_features = np.concatenate(features)
	plot_cluster(X=cat_features[:,0], Y=cat_features[:,1], cluster_center_indices=cci, 
	    figname=args.output+'.tICA_1st_2nd.png', 
	    x_label='tIC 1', y_label='tIC 2', 
	    xmin=floor(cat_features[:,0].min()), xmax=ceil(cat_features[:,0].max()),
	    ymin=floor(cat_features[:,1].min()), ymax=ceil(cat_features[:,1].max()),
	    c_map='winter', cc_color='red')
    elif args.decomposer == 'PCA':
	cat_features = np.concatenate(features)
	plot_cluster(X=cat_features[:,0], Y=cat_features[:,1], cluster_center_indices=cci, 
	    figname=args.output+'.PCA_1st_2nd.png', 
	    x_label='PC 1', y_label='PC 2',
	    xmin=floor(cat_features[:,0].min()), xmax=ceil(cat_features[:,0].max()),
	    ymin=floor(cat_features[:,1].min()), ymax=ceil(cat_features[:,1].max()),
	    c_map='winter', cc_color='red')


if __name__ == '__main__':
    main()
