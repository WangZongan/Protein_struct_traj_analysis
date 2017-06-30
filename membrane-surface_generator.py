#!/usr/bin/env python
'''
Given a membrane thickness, generate the inside and outside surfaces of the membrane. 
Example:
0         1         2         3         4         5 
012345678901234567890123456789012345678901234567890123456789
HETATM 2703  N   DUM  2703     -12.000 -20.000 -14.750                          
HETATM 2703  O   DUM  2703     -12.000 -20.000  14.750                          
HETATM 2704  N   DUM  2704     -12.000 -18.000 -14.750                          
HETATM 2704  O   DUM  2704     -12.000 -18.000  14.750
'''

__author__  = 'Wang Zongan' 
__version__ = '2015-05-18'


def main():
    import sys, string 
    import numpy as np

    thick = sys.argv[1]
    half  = float(thick)/2. 

    f = open('membrane_'+str(thick)+'.pdb','w')
    
    num = 1
    for x in np.arange(-30,31,2):
	for y in np.arange(-30,31,2):
	    print >> f, '%s%5i%3s%6s %5i    %8.3f%8.3f%8.3f' % ('HETATM',num,'N','DUM',num,x,y,-half)
	    print >> f, '%s%5i%3s%6s %5i    %8.3f%8.3f%8.3f' % ('HETATM',num,'O','DUM',num,x,y, half)
	    num += 1

    f.close()

if __name__ == '__main__':
    main() 
