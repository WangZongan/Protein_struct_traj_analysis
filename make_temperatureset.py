#!/usr/bin/env python
import sys
import numpy as np 

Tlow = float(sys.argv[1])
Thigh = float(sys.argv[2])
n_systems = int(sys.argv[3])

assert Thigh > Tlow
assert n_systems > 1

with open('temperature','w') as f:
    print >> f, ','.join(np.linspace(Tlow, Thigh, n_systems).astype('string'))
