import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
# markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
colors = itertools.cycle(('r', 'g')) 
colors2 = itertools.cycle(('b', 'tab:orange')) 

L_cutoffs = list(range(4,10))
thresh_max = 30
for sim_num in range(1, 5):
    L_path = ['/home/agheerae/results/avg_sims/pertnet/', '/home/agheerae/results/sim'+str(sim_num)+'/pertnet/']
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
    ax = []
    p1, p2 = [], []
    for i in range(len(sims_inside)):
        color = next(colors)
        color2 = next(colors2)
        # if i ==0:
        #     lines.append(mlines.Line2D([], [], color=color, label='Average of 4'))
        # else:
        #     lines.append(mlines.Line2D([], [], color=color, label='Simulation '+str(i)))
        for cutoff in L_cutoffs:
            ax.append(f.add_subplot(2, 3, cutoff-3))
            # marker = next(markers)
            ax[cutoff-4].set_title('Cutoff '+str(cutoff)+ 'Å')
            p1.append(ax[cutoff-4].bar(np.asarray(list(range(thresh_max)))+i*0.5, sims_inside[i][cutoff-3, :thresh_max], width=0.5, color=color))
            p2.append(ax[cutoff-4].bar(np.asarray(list(range(thresh_max)))+i*0.5, sims_outside[i][cutoff-3, :thresh_max], width=0.5, color=color2,
                bottom=sims_inside[i][cutoff-3, :thresh_max]))
            plt.xlabel('Threshold')
            plt.ylabel('Number of nodes')
            ax[cutoff-4].set_yscale('log')
            # if i == 0:      
            #     lines.append(mlines.Line2D([], [], color='black', marker=marker, label='Cutoff '+str(cutoff)+'A'))

    f.legend((p1[0], p2[0], p1[-1], p2[-1]), ('Inside the allosteric pathway (average)', 'In the perturbation network (average)', 'Inside the allosteric pathway (sim '+str(sim_num)+')', 'In the perturbation network (sim '+str(sim_num)+')'), loc='upper right', fontsize=8)

    plt.subplots_adjust(hspace=0.5, wspace=0.5, top = 0.8, left=0.1, right=1, bottom=0.1)
    plt.savefig('/home/agheerae/results/bar_sims_nodes_'+str(sim_num)+'_log.pdf')






