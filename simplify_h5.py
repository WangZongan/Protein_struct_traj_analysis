#!/usr/bin/env python
import sys
import numpy as np
import tables as tb 

default_filter = tb.Filters(complib='zlib', complevel=5, fletcher32=True)

def create_array(grp, nm, obj=None):
    return new_t.create_earray(grp, nm, obj=obj, filters=default_filter)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_h5', help='Input simulation file')
    parser.add_argument('output_h5', help='Output similified simulation file')
    args = parser.parse_args()

    t = tb.open_file(args.input_h5) 
    seq = t.root.input.sequence
    pos = t.root.output.pos
    time = t.root.output.time
    potential = t.root.output.potential 

    global new_t 
    new_t = tb.open_file(args.output_h5,'w')

    input = new_t.create_group(new_t.root, 'input')
    new_seq = np.array([s for s in seq])
    create_array(input, 'sequence', obj=new_seq)

    nmodel, nsystem, natom = pos.shape[0:3]
    output = new_t.create_group(new_t.root, 'output')
    new_time = np.zeros(nmodel)
    new_pos = np.zeros((nmodel, nsystem, natom, 3), dtype='f4')
    new_potential = np.zeros((nmodel,nsystem), dtype='f4')
    for i in range(nmodel):
        new_pos[i] = pos[i] 
        new_time[i] = time[i] 
        new_potential[i] = potential[i] 
    create_array(output, 'pos', obj=new_pos)
    create_array(output, 'time', obj=new_time)
    create_array(output, 'potential', obj=new_potential)

    t.close()
    new_t.close()

if __name__ == '__main__':
    main()
