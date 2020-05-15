#!/usr/bin/env python3

import sys,os

def main():
    if len(sys.argv) < 2:
        print("usage : `{}' <csv> [separator=',']".format(sys.argv[0]))
        sys.exit(1)
    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print("ERROR: cannot open '{}'".format(infile))
        sys.exit(1)
    sep = ','
    if len(sys.argv) > 2:
        sep = sys/argv[2]
        if sep == "\\t":
            sep = '\t'
    base = infile.split('.')
    if len(base) > 1:
        outfile = '.'.join(base[:-2]+[base[-2]+"_unique",base[-1]])
        mappingfile = '.'.join(base[:-2]+[base[-2]+"_mapping",base[-1]])
    else:
        outfile=infile+"_unique"
        mappingfile=infile+"_mapping"
    with open(infile, 'rt') as f:
        data = {}
        amb_map = {}
        rev_map = {}
        for idx,line in enumerate(f):
            l=[x.strip() for x in line.split(sep)]
            if idx == 0:
                fields = len(l)
                header = l
                continue
            if len(l) != fields:
                print("ERROR: {} != {}".format(len(l),fields))
            if l[0] not in data:
                data[l[0]] = {}
            if l[1] not in  data[l[0]]:
                data[l[0]][l[1]] = []
            data[l[0]][l[1]].append([x[4:] for x in l[2:]])
    ambig_allele,unambig_allele = 0,0
    ambig_pair,unambig_pair = 0,0
    with open(outfile,'wt') as f:
        for ID in data:
            if len(data[ID]) != 2:
                print("ERROR: should be 2 alleles, got {}".format(len(data[ID])))
            both_unambig = True
            resolved = []
            for allele in data[ID]:
                if len(data[ID][allele]) == 1:
                    unambig_allele+=1
                    resolved.append(data[ID][allele][0])
                else:
                    both_unambig = False
                    ambig_allele+=1
                    r = []
                    for loc_i in range(fields-2):
                        c = []
                        for amb_i in range(len(data[ID][allele])):
                            l = data[ID][allele][amb_i][loc_i]
                            if l != '':
                                c.append(l)
                        if len(c) > 1:
                            s = '|'.join(c)
                            if s not in rev_map:
                                new_id = "{}*ID{}".format(header[loc_i+2][4:],len(amb_map))
                                rev_map[s] = new_id
                                amb_map[new_id] = s 
                                r.append(new_id)
                            else:
                                r.append(rev_map[s])
                        elif len(c) == 0:
                            r.append('')
                        else:
                            r.append(c[0])
                    resolved.append(r)
            if both_unambig:
                unambig_pair+=1
            else:
                ambig_pair+=1
            mux = list(zip(*resolved))
            print(sep.join([ID]+[sep.join(m) for m in mux]), file=f)
    with open(mappingfile,'wt') as f:
        for i in amb_map:
            print("{}={}".format(i,amb_map[i]), file=f)
        
    print("{} ambiguous, {} unambiguous alleles".format(ambig_allele, unambig_allele))
    print("{} ambiguous, {} unambiguous pairs".format(ambig_pair, unambig_pair))



if "__main__" == __name__:
    main()

