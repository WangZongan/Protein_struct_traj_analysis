#!/usr/bin/env python
'''
This program is to perform the clustering analysis on given PDB file. 

Cluster analysis is the grouping of items into clusters based on the similarity of the items to each other.

The following four clustering approaches are implemented in Bio.Cluster:
1. Hierarchical clustering (pairwise centroid-, single-, complete-, and average-linkage); 
2. k-means, k-medians, and k-medoids clustering;
3. Self-Organizing Maps;
4. Principal Component Analysis. 

<Biopython Tutorial and Cookbook, Edi. 2015.12.16 (Biopython 1.66+)>

The data to be clustered are represented by a nxm Numerical Python array data. 
Within the context of gene expression data clustering, typically the rows correspond to different genes
whereas the columns correspond to different experimental conditions. 
The clustering algorithms in Bio.Cluster can be applied both to rows (genes) and to columns (experiments).
'''

__author__  = 'Wang Zongan'
__version__ = '2016-01-19'

import os
import sys
import string
import numpy as np
np.set_printoptions(precision=4)

import cPickle as cp

import Bio.PDB 

backbone = ['C','CA','N'] 
backbone_O = ['C','CA','N','O']  


def calculate_RMSD_matrix(structure, selection): 
    '''
    Return a matrix of RMSD values, in which each RMSD value 
    stands for the difference between 2 models in the structure. 
    '''
    nmodel = len(structure)
    rmsds = np.zeros((nmodel,nmodel))
    def get_aligned_atoms(model_index): 
        atoms = []
        for res in structure[model_index].get_residues():
            if res.get_id()[0] == ' ':  # exclude HET residues and water molecures 
                if selection == 'CA':
                    atoms.append(res['CA'])
                elif selection == 'backbone':
                    for atom_name in backbone:
                        atoms.append(res[atom_name]) 
                elif selection == 'backbone_O': 
                    for atom_name in backbone_O: 
                        atoms.append(res[atom_name]) 
                elif selection == 'all':
                    for atom in res:
                        atoms.append(atom)
        return atoms 
    for i in range(nmodel):
        ref_atoms = get_aligned_atoms(i)
        for j in range(nmodel)[:i]:
            com_atoms = get_aligned_atoms(j) 
            assert len(ref_atoms) == len(com_atoms)
            super_imposer = Bio.PDB.Superimposer() # initiate the superimposer
            super_imposer.set_atoms(ref_atoms,com_atoms)
            super_imposer.apply(com_atoms)
            rmsds[i,j] = super_imposer.rms
            rmsds[j,i] = super_imposer.rms
    return rmsds 


def calculate_TM_score():
    '''
    In bioinformatics, the template modeling score or TM-score is a measure of similarity 
    between two protein structures with different tertiary structures. 

    The TM-score indicates the difference between two structures by a score between (0,1], 
    where 1 indicates a perfect match between two structures.
    < Zhang Y and Skolnick J, Proteins 57(2004)702, 
    Scoring function for automated assessment of protein structure template quality >

    Generally scores below 0.20 corresponds to randomly chosen unrelated proteins 
    whereas structures with a score higher than 0.5 assume roughly the same fold.
    < Zhang Y and Skolnick J, Nucleic Acids Res 33(2005)2302, 
    TM-align: a protein structure alignment algorithm based on the TM-score >

    The TM-score is designed to be independent of protein lengths.

    Formular: 
                       1                                   1
    TM-score = max[ -------- * sum^L_aligned_i (-----------------------) ]
                    L_target                    1 + (di/d0(L_target))^2

    L_target = the length of the target protein
    L_aligned = the length of the aligned region
    d_i = the distance between the ith pair of residues 
    d0(L_target) = 1.24 * (L_target - 15)^(1/3) -1.8

    Note: 'max' exists for example because there are different alignment manners 
    if the target structure and aligned structure have different sequence.  
    '''
    pass 


def calculate_GDT_score():
    pass 


def cluster_centroid_size(cluster_array):
    from collections import Counter
    return Counter(np.array(cluster_array)).most_common()


def hierarchical_clustering_analysis(dist_mat, 
                                    linkage_type='average',  
                                    fcluster_threshold=0.5,
                                    fcluster_criterion='distance'):  
    # fcluster has inconsistent as standard argument for the criterion to choose. 
    # Use distance as argument, to take the cophenetic distance from the linkage matrix Z[:,2]. 
    # You might just use maxclust as criterion if you want to specify the number of clusters. 
    # If you're clustering with single linkage, likely some clusters are singletons (outliers). 

    import scipy.spatial.distance as ssd 
    import scipy.cluster.hierarchy as sch 

    if linkage_type in ['single', 'complete', 'weighted', 'average']:
        condensed_dist_mat = ssd.squareform(dist_mat)
        Z = sch.linkage(condensed_dist_mat, method=linkage_type) 
    elif linkage_type in ['centroid', 'median', 'ward']:
        Z = sch.linkage(dist_mat, method=linkage_type, metric='euclidean')
    
    cluster_array = sch.fcluster(Z,t=dist_mat.max()*fcluster_threshold,criterion=fcluster_criterion)
    return cluster_centroid_size(cluster_array)
    

def pca_analysis():
    pass 


def main():
    import argparse, textwrap
    parser = argparse.ArgumentParser(
        usage ='Use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter, 
        description = textwrap.dedent('''\
            This program employs Bio.Cluster to run clustering analysis 
            on given PDB files containing a trajectory of models.'''))
    parser.add_argument('pdb', help='[required] input PDB file.') 
    parser.add_argument('--measure-type', default='rmsd', type=str,  
        help = textwrap.dedent('''\
            Now, the program only supports RMSD.
            The program will calculate a matrix of scores of this measure type, in which each 
            score stands for the structural difference between two models in the PDB structure.''')) 
    parser.add_argument('--rmsd-selection', default='CA', type=str, 
        help = textwrap.dedent('''\
            This option specifies the selection used to align the models and calculate the RMSD values.
            Available selections are: 
            'CA' (only select CA atoms from the structures), 
            'backbone' (only select backbone atoms: N CA C), 
            'backbone_O' (only select backbone atoms: N CA C O),
            'all' (select all protein atoms). 
            The default value is CA.''')) 

    parser.add_argument('--hierarchical-clustering-analysis', default=False, action='store_true',
        help = 'If turned on, perform Hierarchical Clustering analysis.')
    parser.add_argument('--hierarchical-clustering-analysis-linkage-type', default='average', type=str, 
        help = textwrap.dedent('''\
            The following are methods for calculating the distance between the newly formed cluster u and each v.
            - method='single' assigns 
                    d(u,v)=min(dist(u[i],v[j]))
            for all points i in cluster u and j in cluster v. 
            This is also known as the Nearest Point Algorithm.
            - method='complete' assigns 
                    d(u,v)=max(dist(u[i],v[j]))
            for all points i in cluster u and j in cluster v. 
            This is also known by the Farthest Point Algorithm or Voor Hees Algorithm.
            - method='average' assigns
                    d(u,v)=Sum_{ij}d(u[i],v[j])(|u|*|v|)
            for all points i and j where |u| and |v| are the cardinalities of clusters u and v, respectively. 
            This is also called the UPGMA algorithm. '''))
    parser.add_argument('--hierarchical-clustering-analysis-fcluster-threshold', default=0.5, type=float, 
        help = textwrap.dedent('''\
            The threshold to apply when forming flat clusters. 
            Here, the actual threshold used will be dist_mat.max() * fcluster_threshold. 
            The criterion to use in forming flat clusters is set 'distance' so that the original
            observations in each flat cluster have no greater a cophenetic distance than the threshold. ''')) 

    parser.add_argument('--k-medoids-analysis', default=False, action='store_true', 
        help = 'If turned on, perform K-Medoids analysis.') 
    parser.add_argument('--k-medoids-analysis-nclusters', default=2, type=int,  
        help = 'The number of clusters in K-Medoids analysis. The default value is 2.')

    #parser.add_argument('--pca-analysis', default=False, action='store_true',
    #    help = 'If turned on, perform Principal Component Analysis (PCA) analysis.')

    args = parser.parse_args()

    structure = Bio.PDB.PDBParser(QUIET=True).get_structure('structure',args.pdb)
    nmodel = len(structure) 
    nres = len([res for res in structure[0].get_residues() if res.get_id()[0] == ' '])
    print '=' * 60
    print '%i models detected in the pdb file.' % nmodel
    print '%i residues in the model' % nres
    print '=' * 60

    if args.measure_type == 'rmsd':
        selection = args.rmsd_selection
        distance_matrix = calculate_RMSD_matrix(structure=structure, selection=selection)  

    if args.hierarchical_clustering_analysis:
        distance_matrix_ = distance_matrix
        fc = hierarchical_clustering_analysis(distance_matrix_, 
                        linkage_type=args.hierarchical_clustering_analysis_linkage_type,   
                        fcluster_threshold=args.hierarchical_clustering_analysis_fcluster_threshold)
        print 'hierarchical_clustering_analysis'
        print fc 
        with open(args.pdb+'.hierarchical_clustering_analysis.dat','w') as f:
            for tp in fc: 
                print >> f, '%i %i' % (tp[0],tp[1]) 

    if args.k_medoids_analysis:
        from Bio.Cluster import kmedoids
        '''
        This function implements k-medoids clustering.

        kmedoids(distance, nclusters=2, npass=1, initialid=None)
        
        Arguments:
        - distance: The distance matrix between the elements. There are three ways in which you can pass a distance matrix:
                    1. a 2D Numerical Python array (in which only the left-lower part of the array will be accessed);
                    2. a 1D Numerical Python array containing the distances consecutively;
                    3. a list of rows containing the lower-triangular part of the distance matrix.
                    Examples are:
                    >>> distance = array([[0.0, 1.1, 2.3],
                    ...                   [1.1, 0.0, 4.5],
                    ...                   [2.3, 4.5, 0.0]])
                    (option #1)
                    >>> distance = array([1.1, 2.3, 4.5])
                    (option #2)
                    >>> distance = [array([]),
                    ...             array([1.1]),
                    ...             array([2.3, 4.5])]
                    (option #3)
                    These three correspond to the same distance matrix.
        - nclusters: number of clusters (the 'k' in k-medoids)
        - npass: the number of times the k-medoids clustering algorithm is performed, each time with a different (random) initial condition.
        - initialid: the initial clustering from which the algorithm should start.
                     If initialid is not given, the routine carries out npass repetitions of the EM algorithm, each time starting from a
                     different random initial clustering. If initialid is given, the routine carries out the EM algorithm only once, starting
                     from the initial clustering specified by initialid and without randomizing the order in which items are assigned to
                     clusters (i.e., using the same order as in the data matrix).
                     In that case, the k-means algorithm is fully deterministic.
        Return values:
        - clusterid: array containing the number of the cluster to which each gene/microarray was assigned in the best k-means clustering
                     solution that was found in the npass runs;
        - error: the within-cluster sum of distances for the returned k-means clustering solution;
        - nfound: the number of times this solution was found.

        Returns: clusterid, error, nfound
        '''
        distance_matrix_ = distance_matrix
        clusterid, error, nfound = kmedoids(distance=distance_matrix_, nclusters=args.k_medoids_analysis_nclusters)
        print 'k_medoids_analysis'
        print clusterid
        print cluster_centroid_size(clusterid)
        with open(args.pdb+'.k_medoids_analysis.dat','w') as f:
            for tp in cluster_centroid_size(clusterid):
                print >> f, '%i %i' % (tp[0],tp[1])

if __name__ == '__main__':
    main()

