import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import networkx as nx
marker = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
path1 = '/home/agheerae/results/sim1/nets/threshold/'
path2 = '/home/agheerae/results/sim1/pertnet/'
L_cutoffs = list(range(3,10))
thresh_max = 99

def n_edges(filepath):
    try:
        return len(nx.read_gpickle(filepath).edges())
    except FileNotFoundError:
        return 0

lines = []
f = plt.figure()
for i in L_cutoffs:
    edges_apo, edges_prfar, edges_pertnet = [], [], []
    for j in range(thresh_max+1):
        edges_apo.append(n_edges(path.join(path1, 'apo_'+str(i)+'_'+str(j)+'.p'))) 
        edges_prfar.append(n_edges(path.join(path1, 'prfar_'+str(i)+'_'+str(j)+'.p'))) 
        edges_pertnet.append(n_edges(path.join(path2, 'cutoff_'+str(i), str(i)+'_'+str(j)+'.p')))
    new_marker = next(marker)
    plt.plot(range(thresh_max+1), edges_apo, c='r', marker=new_marker)
    plt.plot(range(thresh_max+1), edges_prfar, c='g', marker=new_marker)
    plt.plot(range(thresh_max+1), edges_pertnet, c='b', marker=new_marker)

    lines.append(mlines.Line2D([], [], color='black', marker=new_marker, label='Cutoff '+str(i)+'Å'))

lines.append(mlines.Line2D([], [], color='r', label='Apo "static"'))
lines.append(mlines.Line2D([], [], color='g', label='PRFAR "static"'))
lines.append(mlines.Line2D([], [], color='b', label='Perturbation network'))

plt.xlabel('Threshold')
plt.ylabel('Number of edges in the networks')
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/number_of_edges.pdf')





