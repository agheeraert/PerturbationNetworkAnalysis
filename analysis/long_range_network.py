from os import listdir
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import matplotlib
import pickle as pkl
import pandas as pd

H = False	
BASE_FOLDER = '/home/agheerae/results/avg_sims/pertnet/'
OUT_FOLDER = '/home/agheerae/results/cutoff_analysis/LRN/15'
tuple_residues = [(['R', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y', 'H'], 'p'),
                (['I', 'L', 'V', 'M', 'F', 'W', 'C', 'P', 'G', 'A'], 'h')]
aa_list = [aa for tup in tuple_residues for aa in tup[0]]
q = len(aa_list)
aa_to_id = dict(zip(aa_list, range(q)))

if __name__ =='__main__':
    for cutoff in range(3, 10):
        folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
        if cutoff == 3 and not H: #FLAG3
            L_threshs = sorted([float(file.split('_')[1].split('.')[0].replace('-','.')) for file in listdir(folder) if file[-1]=='p'])
        else:
            L_threshs = sorted([int(file.split('_')[1].split('.')[0]) for file in listdir(folder) if file[-1]=='p'])
        interactions_count = len(L_threshs)*[0] 
        c_tot = (len(L_threshs))*[0]
        for threshold in L_threshs:
            if cutoff == 3 and not H: #FLAG3
                thresh_str = str(threshold).replace('.', '-')
                if threshold == 1.0:
                    thresh_str = '1'
                elif threshold == 0.0:
                    thresh_str = '0'  
                threshold = int(threshold*10)
            else:
                thresh_str = str(threshold)
            if str(cutoff)+'_'+thresh_str+'.p' in listdir(folder):
                net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+thresh_str+'.p'))
                L_edge_remove = []
                for u, v in net.edges():
                    if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) <= 15 or u[-1] != v[-1]:
                        L_edge_remove.append((u,v))
                        interactions_count[threshold]+=1
                    c_tot[threshold]+=1
                net.remove_edges_from(L_edge_remove)
                net.remove_nodes_from(list(nx.isolates(net)))
            nx.write_gpickle(net, jn(OUT_FOLDER, 'LRN_'+str(cutoff)+'_'+str(threshold)+'.p'))
        interactions_count = [a/b*100 for a,b in zip(interactions_count, c_tot) if b !=0]
        (c_tot, interactions_count)
        lines = []
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.tick_params(labelsize='large')
        ax.set_ylim(0, 101)
        lines.append(mlines.Line2D([], [], color='black', label="Percentage of interactions in the Long Range Network"))
        ax.plot(L_threshs, interactions_count, c='black')
        ax2 = ax.twinx()
        ax2.tick_params(labelsize='large')
        ax2.semilogy(L_threshs, c_tot, color='black', linestyle='dotted')
        ax2.set_ylim(1, 40000)
        lines.append(mlines.Line2D([], [], color='black', label='Total', linestyle=':'))
        ax.grid(linestyle='--')
        #f.legend(handles=lines, loc='upper center', fontsize=8, bbox_to_anchor=(0.5, 1),
        #    ncol=2, fancybox=True, shadow=True)
        plt.tight_layout()
        plt.savefig(jn(OUT_FOLDER, 'long_range_network'+str(cutoff)+'.png'))
        #df = pd.DataFrame(L_top_unknown, index=range(last_cutoff+1))
        #df.to_csv(OUT_FOLDER+'top_unknown_article_2.csv')
