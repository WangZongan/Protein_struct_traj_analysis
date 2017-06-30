#!/usr/bin/env python
import os, sys
import numpy as np


def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
    import argparse
    import re

    if re.match('^([0-9]+(-[0-9]+)?)(,[0-9]+(-[0-9]+)?)*$', s) is None:
        raise argparse.ArgumentTypeError('segments must be of the form 10-30,45,72-76 or similar')

    def parse_seg(x):
        atoms = x.split('-')
        if len(atoms) == 1:
            return np.array([int(atoms[0])])
        elif len(atoms) == 2:
            return np.arange(int(atoms[0]),1+int(atoms[1]))  # inclusive on both ends
        else:
            raise RuntimeError('the impossible happened.  oops.')

    ints = np.concatenate([parse_seg(a) for a in s.split(',')])
    ints = np.array(sorted(set(ints)))   # remove duplicates and sort
    return ints


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'make z flat bottom input file.', 
        usage ='use "%(prog)s --help" for more information')

    parser.add_argument('nres', type=int, help='[required] nres.') 
    parser.add_argument('z0', type=float, help='[required] z0.') 
    parser.add_argument('radius', type=float, help='[required] radius.') 
    parser.add_argument('const', type=float, help='[required] const.') 

    parser.add_argument('--restraint-group', default=[], action='append', type=parse_segments,
        help = 'Multiple residue groups may be specified by giving the --residue-group flag multiple times. ' +
               'Note: residue indices start from 0.' )
    args = parser.parse_args()
    
    if args.restraint_group:
        with open('z_wall.dat', 'w') as f:
            print >> f, 'residue z0 radius spring_constant'
            for rg in args.restraint_group:
                for i in rg:
                    print >> f, '%i %.1f %.1f %.1f' % (i, args.z0, args.radius, args.const)
    else:
        with open('z_wall.dat', 'w') as f:
            print >> f, 'residue z0 radius spring_constant'
            for i in range(nres):
                print >> f, '%i %.1f %.1f %.1f' % (i, args.z0, args.radius, args.const)

    return True


if __name__ == '__main__':
    main()


