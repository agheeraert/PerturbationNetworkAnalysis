import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
marker = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 

path1 = '/home/agheerae/results/all_simu/cutoff/'
path2 = '/home/agheerae/results/sim1/cutoff/'
excel1 = path.join(path1, 'allosteric_path1s_edges.xlsx')
excel2 = path.join(path2, 'allosteric_path1s_edges.xlsx')
L_cutoffs = list(range(3,10))
thresh_max = 40

all_simu = pd.read_excel(excel1).values
sim1 = pd.read_excel(excel2).values

# sim1_Y, all_sim_Y = [], []
lines = []
f = plt.figure()
for i in L_cutoffs:
    new_marker = next(marker)
    plt.plot(range(thresh_max+1), sim1[i-3,:thresh_max+1], c='g', marker=new_marker)
    plt.plot(range(thresh_max+1), all_simu[i-3,:thresh_max+1], c='r', marker=new_marker)
    lines.append(mlines.Line2D([], [], color='black', marker=new_marker, label='Cutoff '+str(i)))

lines.append(mlines.Line2D([], [], color='g', label='Simulation 1'))
lines.append(mlines.Line2D([], [], color='r', label='Average of 4'))

plt.xlabel('Threshold')
plt.ylabel('Number of edges inside the allosteric pathway')
plt.ylim(0, 2000)
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/comparison_sim_edges_tot.pdf')






