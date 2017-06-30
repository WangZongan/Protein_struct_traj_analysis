#!/usr/bin/env python
'''
This program is to perform the clustering analysis on a given sample. 
Cluster analysis is the grouping of items into clusters based on the similarity of the items to each other.

This program takes only a numpy array as input, which is in the shape (n_samples, n_features).

'''

__author__  = 'Wang Zongan'
__version__ = '2016-08-04'

import os
import sys
import string
import numpy as np
import cPickle as cp

from collections import Counter
from sklearn.cluster import KMeans

def cluster_centroid_size(cluster_array):
    return Counter(np.array(cluster_array)).most_common()

def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage ='Use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter, 
        description = textwrap.dedent('''\
            This program employs sklearn.cluster to run clustering analysis 
            on a given numpy array in the shape (n_samples, n_features).'''))
    parser.add_argument('dat', help='[required] input data file (in pkl format).') 
    parser.add_argument('--n-clusters', default=5, type=int,  
        help = textwrap.dedent('''\
                Number of clusters. Default value: 5. ''')) 
    args = parser.parse_args()

    with open(args.dat,'r') as f:
        X = cp.load(f)

    clusterer = KMeans(n_clusters=args.n_clusters)
    clusters  = clusterer.fit(X)
    transformed_X = clusters.transform(X)
    
    # cluster center representative index
    ccri = [np.where(transformed_X[:,i] == transformed_X[:,i].min())[0][0] for i in range(args.n_clusters)]

    with open(args.dat+'.cluster_labels.pkl','w') as f:
        cp.dump(clusters.labels_, f, -1)
    with open(args.dat+'.cluster_center_indices.pkl','w') as f:
        cp.dump(np.array(ccri), f, -1)


if __name__ == '__main__':
    main()

