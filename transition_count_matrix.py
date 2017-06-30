#!/usr/bin/env python
'''
This program is to count the transitions between different states 
and thus to construct the transition count matrix. 

This program takes only an int-type numpy array as input, 
which is in the shape (n_samples,) with each integer refers to a state label.

'''

__author__  = 'Wang Zongan'
__version__ = '2016-08-04'

import os
import sys
import string
import numpy as np
import cPickle as cp

from collections import Counter


def cluster_centroid_size(label_array):
    return Counter(np.array(label_array)).most_common()


def find_transitions(label_array, i, j, n):
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


def main():
    datafile = sys.argv[1]
    lag_time = int(sys.argv[2])

    with open(datafile, 'r') as f:
        labels = cp.load(f)
    
    cluster_sizes = sorted(cluster_centroid_size(labels))
    n_clusters = len(cluster_sizes)

    transition_count_mat = np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        transition_count_mat[i][i] = cluster_sizes[i][1]
        for j in range(i+1,n_clusters):
            transition_count_mat[i][j] = find_transitions(labels,i,j,lag_time)
            transition_count_mat[j][i] = find_transitions(labels,j,i,lag_time)

    print labels
    print transition_count_mat
    with open(datafile+'.transition_count_matrix.pkl','w') as f:
        cp.dump(transition_count_mat, f, -1)


if __name__ == '__main__':
    main()

