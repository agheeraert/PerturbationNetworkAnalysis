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
for is_log in [True, False]:
    for parameter in ['edges', 'nodes']:
        for sim_num in range(1, 5):
            L_path = ['/home/agheerae/results/avg_sims/pertnet/', '/home/agheerae/results/sim'+str(sim_num)+'/pertnet/']
            L_excel = []
            L_excels = []

            for filepath in L_path:
                L_excel.append(path.join(filepath, 'allosteric_path1_'+parameter+'.xlsx'))
                L_excels.append(path.join(filepath, 'allosteric_path1s_'+parameter+'.xlsx'))

            sims = []
            sims_tot = []
            for excel in L_excel:
                sims.append(pd.read_excel(excel).values)
            for excels in L_excels:
                sims_tot.append(pd.read_excel(excels).values)

            sims_inside = []

            for i in range(len(sims)):
                sims_inside.append(sims[i]*sims_tot[i])
            lines = []
            f = plt.figure()
            ax = []
            p1, p2 = [], []
            for i in range(len(sims_inside)):
                color = next(colors)
                color2 = next(colors2)
                if i ==0:
                    lines.append(mlines.Line2D([], [], color=color, label='Inside the allosteric pathway (Average of 4)'))
                    lines.append(mlines.Line2D([], [], color=color2, label='In the perturbation network (Average of 4)'))
                else:
                    lines.append(mlines.Line2D([], [], color=color, label='Inside the allosteric pathway (Simulation '+str(sim_num)+')'))
                    lines.append(mlines.Line2D([], [], color=color2, label='In the perturbation network (Simulation '+str(sim_num)+')'))
                for cutoff in L_cutoffs:
                    ax.append(f.add_subplot(2, 3, cutoff-3))
                    # marker = next(markers)
                    ax[cutoff-4].set_title('Cutoff '+str(cutoff)+ 'Å')
                    p1.append(ax[cutoff-4].plot(range(thresh_max), sims_inside[i][cutoff-3,:thresh_max], c=color))
                    p2.append(ax[cutoff-4].plot(range(thresh_max), sims_tot[i][cutoff-3,:thresh_max], c=color2))
                    plt.xlabel('Threshold')
                    plt.ylabel('Number of '+parameter)
                    if is_log:
                        ax[cutoff-4].set_yscale('log')

            f.legend(handles=lines, loc='upper left', fontsize=8)
            plt.subplots_adjust(hspace=0.5, wspace=0.5, top = 0.8, left=0.1, right=1, bottom=0.1)
            if is_log:
                plt.savefig('/home/agheerae/results/ratio_sims_'+parameter+'_'+str(sim_num)+'_log.pdf')
            else:
                plt.savefig('/home/agheerae/results/ratio_sims_'+parameter+'_'+str(sim_num)+'.pdf')






