#!/usr/bin/env python
__author__  = 'Wang Zongan'
__version__ = '2016-06-24'

import os
import sys 
import string
import numpy as np
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def plot_contour(X, Y, xname, yname, outputfilename, plot_format='png'):
    fig = figure(figsize=(10,10))
    title('Contour map: %s --> %s' % (xname, yname),fontsize=25)
    plot(X, Y, linestyle='--', marker='o', color='b', markersize=1)
    tick_params(axis='both', which='major', labelsize=15)
    legend(loc='upper left', fontsize='x-large')
    xlabel(xname, fontsize=15)
    ylabel(yname, fontsize=15)
    lmt = (int(np.maximum(X.max(),Y.max())/10.)+1)*10
    xlim(0,lmt)
    ylim(0,lmt)
    grid()
    tight_layout()
    savefig(outputfilename, format=plot_format)


def main():
    X = np.loadtxt(sys.argv[1])
    Y = np.loadtxt(sys.argv[2])
    xname = sys.argv[3]
    yname = sys.argv[4]
    outputfilename = sys.argv[5] 
    plot_contour(X, Y, xname, yname, outputfilename)
        
if __name__ == '__main__':
    main()
