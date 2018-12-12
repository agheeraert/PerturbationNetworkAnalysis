import pandas as pd
from os import path
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
# markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
colors = itertools.cycle(('r', 'g', 'b','olive', 'purple', 'navy', 'brown')) 
colors2 = itertools.cycle(('c', 'y', 'm', 'tab:orange', 'pink', 'lime', 'turquoise')) 

L_cutoffs = list(range(4,10))
thresh_max = 30
for is_log in [True, False]:
    for parameter in ['edges', 'nodes']:
        # for time_win in range(0, 6):
        #     L_path = ['/home/agheerae/results/sim1/pertnet/', '/home/agheerae/results/timewin_'+str(time_win*100)+'/pertnet/']
        L_path = ['/home/agheerae/results/sim1/pertnet/']
        for time_win in range(0, 6):
            L_path.append('/home/agheerae/results/timewin_'+str(time_win*100)+'/pertnet/')
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
        ax1, ax2 = [], []
        for i in range(len(sims_inside)):
            color = next(colors)
            color2 = next(colors2)
            if i ==0:
                lines.append(mlines.Line2D([], [], color=color, label='Inside the allosteric pathway (0-100ns)'))
                lines.append(mlines.Line2D([], [], color=color2, label='In the perturbation network (0-100ns)'))
            else:
                lines.append(mlines.Line2D([], [], color=color, label='Inside the allosteric pathway ('+str((i-1)*10)+'-'+str((4+i)*10)+'ns)'))
                lines.append(mlines.Line2D([], [], color=color2, label='In the perturbation network ('+str((i-1)*10)+'-'+str((4+i)*10)+'ns)'))
            for cutoff in L_cutoffs:
                ax1.append(f.add_subplot(2, 3, cutoff-3))
                # marker = next(markers)
                ax1[cutoff-4].set_title('Cutoff '+str(cutoff)+ 'Å')
                ax1[cutoff-4].plot(range(thresh_max), sims_inside[i][cutoff-3,:thresh_max], c=color)
                ax1[cutoff-4].set_xlabel('Threshold')
                ax1[cutoff-4].set_ylabel(ylabel='Number of '+parameter)
                ax1[cutoff-4].tick_params(axis='y', labelcolor=color)
                ax2.append(ax1[-1].twinx())
                ax2[cutoff-4].plot(range(thresh_max), sims_tot[i][cutoff-3,:thresh_max], c=color2)
                ax2[cutoff-4].tick_params(axis='y', labelcolor=color2)
                if is_log:
                    ax1[cutoff-4].set_yscale('log')
                    ax2[cutoff-4].set_yscale('log')

        f.legend(handles=lines, loc='upper left', fontsize=8)
        plt.subplots_adjust(hspace=0.5, wspace=0.75, top = 0.8, left=0.1, right=1, bottom=0.1)
        if is_log:
            plt.savefig('/home/agheerae/results/diffscale_'+parameter+'_timewin_'+'log.pdf')
        else:
            plt.savefig('/home/agheerae/results/diffscale_'+parameter+'_timewin'+'.pdf')






