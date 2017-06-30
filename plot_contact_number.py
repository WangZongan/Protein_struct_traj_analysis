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


def plot_contact_number_evolution(
        contacts, native_contacts, pdbfile, contact_type, cutoff, plot_format='png'):

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


def plot_heat_map_contact_number(
        contacts, native_contacts, pdbfile, contact_type, cutoff,
        plot_color_map='jet', plot_format='png'):  

    heatmap, xedges, yedges = np.histogram2d(
            contacts, native_contacts, 
            bins=[np.linspace(contacts.min(), contacts.max(), 51),
                  np.linspace(native_contacts.min(), native_contacts.max(), 51)],
            normed=True)
    
    extent = (contacts.min(), contacts.max(), native_contacts.min(), native_contacts.max())

    fig = figure(figsize=(12,10))
    suptitle('%s Contact Map: cutoff=%s'%(contact_type, str(cutoff)),fontsize=25,x=0.435,y=0.98)
    title(pdbfile, fontsize=20)
    my_cmap = get_cmap(plot_color_map)
    my_cmap.set_under('w')

    imshow(heatmap, origin='low', cmap=my_cmap, extent=extent, aspect='auto')
    xlim(contacts.min(), contacts.max())
    ylim(native_contacts.min(), native_contacts.max())

    tick_params(axis='both', which='major', labelsize=15)
    fig.text(0.4, 0.03, '# total contacts', ha='center', va='center', fontsize=25)
    fig.text(0.02, 0.5, '% native contacts', ha='center', va='center', fontsize=25, rotation='vertical')
    grid()

    cb = colorbar()
    #cb.set_label()
    tight_layout(rect=[0.03, 0.05, 1, 0.95]) # default is (0, 0, 1, 1) [left, bottom, right, top]
    savefig('%s.heat_map_contact_number.%s.%s.%s' % 
            (pdbfile, contact_type, str(cutoff), plot_format), 
        format=plot_format)  


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Calculate and plot contact map for a given pdb file, ' +
                      'which can contain a trajectory of structures.',
        usage ='use "%(prog)s --help" for more information')
    parser.add_argument('dat', help='[required] input DAT file')
    parser.add_argument('--contact-type', default='CA', type=str,
        help = 'Contact type, options are CA, CB, CM, which stand for CA-CA contact, '+
               'CB-CB contact, and sidechain CM-CM contact, respectively. '+
               'Now all the contact types supported are distance cutoff contacts.' +
               'The default value is CA.')
    parser.add_argument('--contact-cutoff', default=7.5, type=float,
        help = 'Suggested values for CA, CB, and CM contacts are 7.5, 7.0, 6.5 A, respectively. ' +
               'The default value is 7.5 for the default CA contact.')
    parser.add_argument('--plot-color-map', default='jet', type=str, 
        help = "Color map used in plotting, jet by default. " +
               "Any color map supported by matplotlib can be used. " + 
               "Examples are: 'Blues', 'GnBu', 'BrBG', 'gist_rainbow', etc. " + 
               "(Ref: http://matplotlib.org/xkcd/examples/color/colormaps_reference.html)")
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
