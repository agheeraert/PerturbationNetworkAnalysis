import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
colors = itertools.cycle(('r', 'g', 'b', 'tab:orange', 'tab:purple')) 
colors2 = itertools.cycle(('y', 'm', 'c', 'tab:pink', 'tab:olive')) 

L_cutoffs = list(range(3,10))
thresh_max = 30

L_path = ['/home/agheerae/results/avg_sims/pertnet/', '/home/agheerae/results/sim1/pertnet/', '/home/agheerae/results/sim2/pertnet/', 
'/home/agheerae/results/sim3/pertnet/', '/home/agheerae/results/sim4/pertnet/']
L_excel = []
L_excels = []

for filepath in L_path:
    L_excel.append(path.join(filepath, 'allosteric_path1_nodes.xlsx'))
    L_excels.append(path.join(filepath, 'allosteric_path1s_nodes.xlsx'))

sims = []
sims_tot = []
for excel in L_excel:
    sims.append(pd.read_excel(excel).values)
for excels in L_excels:
    sims_tot.append(pd.read_excel(excels).values)

sims_inside = []
sims_outside = []

for i in range(len(sims)):
    sims_inside.append(sims[i]*sims_tot[i])
    sims_outside.append((1-sims[i])*sims_tot[i])
lines = []
f = plt.figure()
ax = f.add_subplot(111)
for i in range(len(sims_inside)):
    color = next(colors)
    color2 = next(colors2)
    # if i ==0:
    #     lines.append(mlines.Line2D([], [], color=color, label='Average of 4'))
    # else:
    #     lines.append(mlines.Line2D([], [], color=color, label='Simulation '+str(i)))
    for cutoff in L_cutoffs:
        marker = next(markers)
        p1 = ax.bar(np.asarray(list(range(thresh_max)))+i*0.2, sims_inside[i][cutoff-3, :thresh_max], width=0.2, color=color)
        p2 = ax.bar(np.asarray(list(range(thresh_max)))+i*0.2, sims_outside[i][cutoff-3, :thresh_max], width=0.2, color=color2,
             bottom=sims_inside[i][cutoff-3, :thresh_max])
        # if i == 0:      
        #     lines.append(mlines.Line2D([], [], color='black', marker=marker, label='Cutoff '+str(cutoff)+'A'))

plt.xlabel('Threshold')
plt.ylabel('Total number of nodes in the perturbation network')
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/bar_sims_nodes.pdf')






