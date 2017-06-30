#!/usr/bin/env python
'''
Given pdb trajectory, perform FAST analysis.. 
'''
__author__  = 'Wang Zongan'
__version__ = '2016-08-29'

import os
import sys
import string
import numpy as np
import cPickle as cp
import mdtraj as md

import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def calculate_rmsd_mat_mdtraj(com_traj, ref_traj):
    from msmbuilder.featurizer import RMSDFeaturizer
    return RMSDFeaturizer(ref_traj).partial_transform(com_traj)


def cluster_centroid_size(labels):
    from collections import Counter
    return Counter(np.array(labels)).most_common()


def cluster(X, n_clusters=5):
    '''
    X: (n_samples, n_features)
    '''
    from sklearn.cluster import KMeans
    clusterer = KMeans(n_clusters=n_clusters)
    clusters  = clusterer.fit(X)
    transformed_X = clusters.transform(X) # (n_samples, n_clusters)

    # cluster center index. 
    # The centers returned by KMeans are not necessarily in the samples, 
    # I just let those in the samples closest to the KMeans centers to be the centers. 
    cci = [np.where(transformed_X[:,i] == transformed_X[:,i].min())[0][0] for i in range(n_clusters)]

    return clusters.labels_, np.array(cci) 


def plot_cluster(X, cluster_center_indices, figname, ini_tar_RMSD=3):
    fig = figure(figsize=(10,10))
    hexbin(X[:,0], X[:,1], bins='log', mincnt=1, cmap="Blues")
    scatter(X[cluster_center_indices, 0], X[cluster_center_indices, 1], s=100, c='red')
    xlabel('RMSD to initial structure', fontsize=20)
    ylabel('RMSD to target structure', fontsize=20)
    xlim(0, ini_tar_RMSD) # unit: nm
    ylim(0, ini_tar_RMSD)
    tick_params(axis='both', which='major', labelsize=20)
    grid()
    tight_layout()
    savefig(figname, format='png')
    

def find_transitions(label_array, i, j, n=1):
    '''
    Find the number of transitions from state i to state j in label_array.
    n: lag number
    '''
    assert i in label_array
    assert j in label_array
    i_idx = np.where(label_array == i)[0]
    i_idx_n = i_idx + n
    i_idx_n = i_idx_n[np.where(i_idx_n < len(label_array))]
    i_to_j_idx = i_idx[np.where(label_array[i_idx_n] == j)]
    return len(i_to_j_idx)


def calc_transition_count_mat(labels, n_clusters):
    cluster_sizes = sorted(cluster_centroid_size(labels))
    print 'cluster_sizes: ', cluster_sizes

    transition_count_mat = np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        transition_count_mat[i][i] = cluster_sizes[i][1]
        for j in range(i+1, n_clusters):
            transition_count_mat[i][j] = find_transitions(labels,i,j,1)
            transition_count_mat[j][i] = find_transitions(labels,j,i,1)
    return transition_count_mat


def calc_reward_score_pot(pots, cluster_center_indices, transition_count_mat, alpha=1, n_simulations=30):
    '''
    pots: (n_samples,)
    transition_count_mat: (n_clusters, n_clusters)
    '''
    n_clusters = len(cluster_center_indices)

    pot_max = pots[cluster_center_indices].max()
    pot_min = pots[cluster_center_indices].min()
    pot_max_min = pot_max - pot_min 
    
    c_max = np.diag(transition_count_mat).max()
    c_min = np.diag(transition_count_mat).min()
    c_max_min = c_max - c_min
    c = np.sum(transition_count_mat, axis=0) - np.diag(transition_count_mat) # total transition from state i to any other state

    rewards = np.array([ 
        ( (pot_max - pots[cluster_center_indices[i]]) / pot_max_min + alpha * (c_max - c[i]) / c_max_min ) 
        for i in range(n_clusters) ])
    
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

    print '=' * 50
    print '%6s %8s %10s %16s %8s %6s' % ('idx', '#cluster', 'potential', '#To other states', 'reward', 'sim')
    for i in range(n_clusters):
        print '%6i %8i %10.3f %16i %8.3f %6i' % (
                cluster_center_indices[i], 
                transition_count_mat[i][i],
                pots[cluster_center_indices[i]], 
                c[i],
                rewards[i], 
                sims[i])

    return rewards, sims, c


def main():
    import argparse
    parser = argparse.ArgumentParser(
        usage ='use "python %(prog)s --help" for more information',
        description = 'Featurize the given pdb trajectory into a vectorizable space. ')
    parser.add_argument('pdb', help='[required] input PDB file.')
    parser.add_argument('potential', type=str, help = "[required] The potential file corresponding to the pdb file.")
    parser.add_argument('--initial-pdb', type=str, help = "The initial pdb structure. " + 
        "Caution: The initial pdb should have the same number of atoms in structure with that in pdb.") 
    parser.add_argument('--target-pdb', type=str, help = "The target pdb structure. " +
        "Caution: The target pdb should have the same number of atoms in structure with that in pdb.") 
    parser.add_argument('--n-clusters', default=5, type=int, help = "Number of clusters. Default value: 5.")
    parser.add_argument('--n-simulations', default=30, type=int, 
        help = "Number of parallel simulations in each round of FAST algorithm. Default value: 30. ")
    parser.add_argument('--alpha', default=1., type=float, help = "Number of clusters. Default value: 1.0.")
    parser.add_argument('--output', type=str, default='output', 
        help = 'Output the pkl file of the shape (n_samples, n_features).')
    args = parser.parse_args()

    ## trajectory
    traj = md.load(args.pdb) ; print 'Traj loaded.'
    traj_rmsd_mat = calculate_rmsd_mat_mdtraj(traj, traj)     # (n_samples, n_samples), used as cluster_mat

    ini_traj = md.load(args.initial_pdb); print 'Initial traj loaded.'
    tar_traj = md.load(args.target_pdb) ; print 'Target traj loaded.'
    ini_rmsd_mat  = calculate_rmsd_mat_mdtraj(traj, ini_traj) # (n_samples,)
    tar_rmsd_mat  = calculate_rmsd_mat_mdtraj(traj, tar_traj) # (n_samples,)
    ini_tar_rmsd_mat = np.column_stack((ini_rmsd_mat, tar_rmsd_mat))    # (n_samples, 2)

    ## potential 
    pots = np.loadtxt(args.potential)  # (n_samples, )
    assert len(pots) == len(traj_rmsd_mat)

    ## clustering
    labels, cluster_center_indices = cluster(traj_rmsd_mat, args.n_clusters)
    print 'Labels: ', labels
    print 'Cluster centers: ', cluster_center_indices

    transition_count_mat = calc_transition_count_mat(labels, args.n_clusters)
    print 'Transition count matrix: \n', transition_count_mat

    ## plot
    plot_cluster(ini_tar_rmsd_mat, cluster_center_indices, args.output+'.png')

    ## reward score
    rewards, sims, c = calc_reward_score_pot(
                        pots=pots, 
                        cluster_center_indices=cluster_center_indices, 
                        transition_count_mat=transition_count_mat, 
                        alpha=args.alpha, 
                        n_simulations=args.n_simulations)

    # output dat
    with open(args.output+'.iniRMSD.dat', 'w') as f:
        for i in range(len(ini_rmsd_mat)):
            print >> f, '%8.3f' % (ini_rmsd_mat[i])
    with open(args.output+'.tarRMSD.dat', 'w') as f:
        for i in range(len(tar_rmsd_mat)):
            print  >> f, '%8.3f' % (tar_rmsd_mat[i])

    with open(args.output+'.dat', 'w') as f:
        print >> f, '%6s %8s %10s %12s %16s %8s %6s' % (
                    'idx', '#cluster', 'potential', 'RMSD_target', '#To other states', 'reward', 'sim')
        for i in range(args.n_clusters):
            print >> f, '%6i %8i %10.3f %12.3f %16i %8.3f %6i' % (
                cluster_center_indices[i], 
                transition_count_mat[i][i],
                pots[cluster_center_indices[i]],
                tar_rmsd_mat[cluster_center_indices[i]], 
                c[i],
                rewards[i], 
                sims[i])

    with open(args.output+'.CenterIdx_ClusterSize.dat', 'w') as f:
        for i in range(args.n_clusters):
            print >> f, '%6i %6i' % (cluster_center_indices[i], sims[i])


if __name__ == '__main__':
    main()
