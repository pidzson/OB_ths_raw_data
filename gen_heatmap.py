#!/usr/bin/env pyton3
import sys,os
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("usage: {} <interaction_graph.txt>".format(sys.argv[0]))
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print("ERROR: cannot open '{}'".format(sys.argv[1]))
        sys.exit(1)
    
    with open(sys.argv[1], 'rt') as f:
        m={}
        s=set()
        for idx,line in enumerate(f):
            if idx == 0:
                #skip header
                continue
            l = [x.strip() for x in line.split(' ') if x.strip() != '']
            k = (l[0],l[1])
            if l[6] != 'NA' and float(l[6]) < 0.05:
                m[k] = np.log(float(l[6]))
            elif l[7] != 'NA' and float(l[7]) < 0.05:
                m[k] = -np.log(float(l[7]))
            else:
                m[k] = 0.0
            s.add(k[0])
            s.add(k[1])
    l=list(s)
    l.sort()
    print('X\t'+'\t'.join(l))
    for i in l:
        print(i, end='')
        for j in l:
            if (i,j) in m:
                print('\t{}'.format(m[(i,j)]), end='')
            elif (j,i) in m:
                print('\t{}'.format(m[(j,i)]), end='')
            else:
                print('\t0.0', end='')
        print()


if "__main__" == __name__:
    main()
