#!/usr/bin/env python3
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

def main():
    if len(sys.argv) < 2:
        print("usage: {} <interaction.txt>".format(sys.argv[0]))
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print("ERROR: cannot open '{}'".format(sys.argv[1]))
        sys.exit(1)
    
    G=nx.Graph()
    ints = pd.read_csv(sys.argv[1], delim_whitespace=True, header=0)
    for r in ints.iterrows():
        c='k'
        w=0
        if r[1]['P7'] < 0.05:
            # difference
            c='r'
            #w=np.log(r[1]['OR7'])
            #if w > 0:
            #    print('warning {} {} should differ'.format(r[1]['ID1'], r[1]['ID2']))
            w=np.log(r[1]['P7'])
        elif r[1]['P8'] < 0.05:
            # combined
            c='g'
            #w=np.log(r[1]['OR8'])
            #if w < 0:
            #    print('warning {} {} should combine'.format(r[1]['ID1'], r[1]['ID2']))
            w=-np.log(r[1]['P8'])
        G.add_edge(r[1]['ID1'], r[1]['ID2'], color=c, weight=w)
    pos={}
    for n in G.nodes():
        pos[n] = (0,0)
    pos['DRB1*07:01:01']=(-1,.5)
    pos['DQB1*02:02:01']=(-1,-.5)
    fixed=['DRB1*07:01:01', 'DQB1*02:02:01']
    pos=nx.spring_layout(G, pos=None, fixed=None, iterations=50, k=.1)
    edges = G.edges()
    colors = [G[u][v]['color'] for u,v in edges]
    weights = [G[u][v]['weight']*.2 for u,v in edges]
    weights = [.1 if x == 0 else x for x in weights]
    nx.draw_networkx_labels(G, pos)
    nx.draw(G, pos, edges=edges, edge_color=colors, width=weights)
    plt.draw()
    plt.show()

if "__main__" == __name__:
    main()

