import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
colors = itertools.cycle(('r', 'g', 'b', 'c', 'purple')) 

L_cutoffs = list(range(3,10))
thresh_max = 30

L_path = ['/home/agheerae/results/avg_sims/pertnet/', '/home/agheerae/results/sim1/pertnet/', '/home/agheerae/results/sim2/pertnet/', 
'/home/agheerae/results/sim3/pertnet/', '/home/agheerae/results/sim4/pertnet/']
L_excel = []

for filepath in L_path:
    L_excel.append(path.join(filepath, 'allosteric_path1s_nodes.xlsx'))

sims = []
for excel in L_excel:
    sims.append(pd.read_excel(excel).values)

# sim1_Y, all_sim_Y = [], []
lines = []
f = plt.figure()
for i, sim in enumerate(sims):
    color = next(colors)
    if i ==0:
        lines.append(mlines.Line2D([], [], color=color, label='Average of 4'))
    else:
        lines.append(mlines.Line2D([], [], color=color, label='Simulation '+str(i)))
    for cutoff in L_cutoffs:
        marker = next(markers)
        plt.plot(range(thresh_max), sim[cutoff-3,:thresh_max], c=color, marker=marker) 
        if i == 0:      
            lines.append(mlines.Line2D([], [], color='black', marker=marker, label='Cutoff '+str(cutoff)+'A'))

plt.xlabel('Threshold')
plt.ylabel('Total number of nodes in the perturbation network')
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/all_sims_nodes_tot.pdf')






