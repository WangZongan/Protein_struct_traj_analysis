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
    parser.add_argument('input_h5', help='First input similified simulation file') 
    parser.add_argument('--add-input-h5', default=None, 
        help='Second input similified simulation file')   
    parser.add_argument('output_h5', help='Output similified simulation file')
    args = parser.parse_args() 

    t1 = tb.open_file(args.input_h5) 
    seq1 = t1.root.input.sequence
    pos1 = t1.root.output.pos
    time1 = t1.root.output.time
    potential1 = t1.root.output.potential 

    nres = seq1.shape[0]
    nmodel1, nsystem, natom = pos1.shape[0:3]

    global new_t 
    new_t = tb.open_file(args.output_h5,'w')

    input = new_t.create_group(new_t.root, 'input')
    new_seq = np.array([s for s in seq1])
    create_array(input, 'sequence', obj=new_seq)

    output = new_t.create_group(new_t.root, 'output')

    if args.add_input_h5 is not None:
        t2 = tb.open_file(args.add_input_h5) 
        seq2 = t2.root.input.sequence
        pos2 = t2.root.output.pos
        time2 = t2.root.output.time
        potential2 = t2.root.output.potential 
        
        for i in range(nres):
            assert seq1[i] == seq2[i]
        assert pos1.shape[1:3] == pos2.shape[1:3]  # (nmodel, nsystem, natom, 3)
        nmodel2, nsystem, natom = pos2.shape[0:3]
    
        new_pos = np.zeros((nmodel1+nmodel2, nsystem, natom, 3), dtype='f4')
        new_pos[0:nmodel1] = pos1
        new_pos[nmodel1:]  = pos2 

        new_time = np.zeros(nmodel1+nmodel2)
        new_time[0:nmodel1] = time1
        new_time[nmodel1:]  = time2 + time1[-1]

        new_potential = np.zeros((nmodel1+nmodel2,nsystem), dtype='f4')
        new_potential[0:nmodel1] = potential1
        new_potential[nmodel1:]  = potential2 
        t2.close()
    
    else:
        new_time = np.zeros(nmodel1)
        new_potential = np.zeros((nmodel1,nsystem), dtype='f4')
        new_pos = np.zeros((nmodel1, nsystem, natom, 3), dtype='f4')
        for i in range(nmodel1):
            new_pos[i] = pos1[i] 
            new_time[i] = time1[i] 
            new_potential[i] = potential1[i] 

    create_array(output, 'pos', obj=new_pos)
    create_array(output, 'time', obj=new_time)
    create_array(output, 'potential', obj=new_potential)
  
    t1.close()    
    new_t.close()

if __name__ == '__main__':
    main()
