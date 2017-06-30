#!/usr/bin/env python
__author__  = 'Wang Zongan'
__version__ = '2016-08-10'

import os
import sys 
import string
import numpy as np
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def plot_scatter(X, Y, xname, yname, outputfilename, plot_format='png'):
    fig = figure(figsize=(10,10))
    title('Contour map: %s --> %s' % (xname, yname),fontsize=25)
    hexbin(X, Y, bins='log', mincnt=1, cmap="bone_r")
    scatter(clusterer.cluster_centers_[msm.state_labels_, 0],
                clusterer.cluster_centers_[msm.state_labels_, 1],
                s=1e4 * msm.populations_,       # size by population
                c=msm.left_eigenvectors_[:, 1], # color by eigenvector
                cmap="coolwarm") 
    colorbar(label='First dynamical eigenvector')
    xlabel(xname)
    ylabel(yname)
    tight_layout()
    tick_params(axis='both', which='major', labelsize=15)
    savefig(outputfilename, format=plot_format)


def main():
    X = np.loadtxt(sys.argv[1])
    Y = np.loadtxt(sys.argv[2])
    xname = sys.argv[3]
    yname = sys.argv[4]
    outputfilename = sys.argv[5] 
    plot_scatter(X, Y, xname, yname, outputfilename)
        
if __name__ == '__main__':
    main()
