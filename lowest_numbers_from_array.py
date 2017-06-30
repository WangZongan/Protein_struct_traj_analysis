#!/usr/bin/env python
import os
import sys
import string
import numpy as np


def get_lowest_percent_idx(numarr, percent):
    order = numarr.argsort().argsort()
    tot = len(numarr)
    return np.where(order <= tot*percent/100.)[0]


def main():
    numarr = np.loadtxt(sys.argv[1])
    percent = float(sys.argv[2])
    lowest = get_lowest_percent_idx(numarr, percent)

    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0] + '.lowest%.1fpercent.dat' % percent, 'w') as f:
        for e in lowest.astype(str):
            print >> f, e


if __name__ == '__main__':
    main()



