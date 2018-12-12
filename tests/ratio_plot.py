import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
colors = itertools.cycle(('r', 'g', 'b', 'c', 'y', 'olive', 'm')) 

L_cutoffs = list(range(3,10))
thresh_max = 30

L_path = ['/home/agheerae/results/sim1/pertnet/', '/home/agheerae/results/timewin_0/pertnet/', '/home/agheerae/results/timewin_100/pertnet/', 
'/home/agheerae/results/timewin_200/pertnet/', '/home/agheerae/results/timewin_300/pertnet/', '/home/agheerae/results/timewin_400/pertnet/', '/home/agheerae/results/timewin_500/pertnet/']
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
        lines.append(mlines.Line2D([], [], color=color, label='Simulation 1'))
    else:
        lines.append(mlines.Line2D([], [], color=color, label='Time window '+str((i-1)*100)+'-'+str((i+4)*100)+'ns'))
    for cutoff in L_cutoffs:
        marker = next(markers)
        plt.plot(range(thresh_max), sim[cutoff-3,:thresh_max], c=color, marker=marker) 
        if i == 0:      
            lines.append(mlines.Line2D([], [], color='black', marker=marker, label='Cutoff '+str(cutoff)+'Å'))

plt.xlabel('Threshold')
plt.ylabel('Total number of nodes in the perturbation network')
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/time_windows_nodes_tot.pdf')






