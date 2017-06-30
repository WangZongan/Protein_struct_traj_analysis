#!/usr/bin/env python
import sys
import numpy as np

def main():

    lT=float(sys.argv[1])
    hT=float(sys.argv[2])
    ns=int(sys.argv[3])

    Tlist = np.zeros(ns)
    
    for i in range(ns):
        Tlist[i] = ((1-float(i)/float(ns-1))*np.sqrt(lT)+(float(i)/float(ns-1))*np.sqrt(hT))**2
    print Tlist
    
    with open('templist_%.1f_%.1f_%i.dat'%(lT, hT, ns),'w') as f:
        outputstring = ''
        for t in Tlist:
            outputstring += '%.3f,'%t
        outputstring = outputstring.strip(',')
        f.write(outputstring)

    with open('swapset_%.1f_%.1f_%i.dat'%(lT, hT, ns),'w') as f:
        even = ''
        odd = '' 
        n = 0 
        while 2*n < ns-1:
            even += '%i-%i,' % (2*n, 2*n+1)
            n += 1
        n = 0
        while 2*n+1 < ns-1: 
            odd += '%i-%i,' % (2*n+1, 2*n+2)
            n += 1
        even = even.strip(',')
        odd = odd.strip(',')

        outputstring = '--swap-set %s --swap-set %s' % (even, odd)
        f.write(outputstring)

    return Tlist

if __name__ == '__main__':
    main()
