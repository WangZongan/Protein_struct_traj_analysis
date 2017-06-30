#!/usr/bin/env python
__author__  = 'Wang Zongan'
__version__ = '2016-05-03'

import os
import sys 
import string
import numpy as np
import seaborn as sns
sns.set_style(style='white')
from matplotlib.pyplot import *


def plot_lines(data, outputfilename, 
        plot_title, 
        plot_format='png'):

    fig = figure(figsize=(12,10))
    title('# Contacts V.S. Frame',fontsize=25)
    plot(contacts, color='green', label='# total contacts')
    plot(native_contacts, color='red', label='# native contacts')
    tick_params(axis='both', which='major', labelsize=15)
    grid()
    legend(loc='upper left', fontsize='x-large')
    tight_layout()
    savefig('%s.contact_number_evolutoin.%s.%s.%s' % 
                (pdbfile, contact_type, str(cutoff), plot_format), 
            format=plot_format)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Calculate and plot contact map for a given pdb file, ' +
                      'which can contain a trajectory of structures.',
        usage ='use "%(prog)s --help" for more information')
    parser.add_argument('dat', help='[required] input DAT file')
    parser.add_argument('--plot-format', default='png', type=str,
        help = 'Format of output plot, PNG format by default. ' +
               'Any format supported by matplotlib can be used.')
    args = parser.parse_args()

    dat = np.loadtxt(args.dat)
    contacts, native_contacts = dat[:,0], dat[:,1] 

    plot_contact_number_evolution(
                contacts, 
                native_contacts,
                os.path.basename(args.dat), 
                contact_type=args.contact_type, 
                cutoff=args.contact_cutoff,
                plot_format=args.plot_format)
    
    plot_heat_map_contact_number(
                contacts, 
                native_contacts,
                os.path.basename(args.dat),
                contact_type=args.contact_type,
                cutoff=args.contact_cutoff,
                plot_color_map=args.plot_color_map,
                plot_format=args.plot_format)

        
if __name__ == '__main__':
    main()
