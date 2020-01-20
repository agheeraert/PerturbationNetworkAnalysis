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
from Bio.PDB.Polypeptide import aa3, aa1

BASE_FOLDER = '/home/agheerae/results/avg_sims/nets/threshold/'
OUT_FOLDER = '/home/agheerae/results/cutoff_analysis/w_and_k/by_aa_size/'
thresh = 5
colors = itertools.cycle(('orange', 'g', 'b', 'magenta', 'r', 'black', 'navy'))   # 'purple',  between olive and navy
style = itertools.cycle(('-', '--', ':'))
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))

size_to_aa = {5: ['G'],
           6: ['A'],
           7: ['C', 'S']
           8: ['P', 'T'; 'V']
           9: ['D', 'I', 'L', 'M', 'N']
           10: ['E', 'K', 'Q']
           11: ['H'],
           12: ['F', 'R']
           13: ['Y']
           15: ['W']}

if __name__ =='__main__':
    WK = [] 
    for num, folder in enumerate(['apo', 'prfar']):
        WK.append(np.zeros((7, q)))
        for cutoff in range(3, 10):
            subfolder = jn(BASE_FOLDER, folder, 'cutoff_'+str(cutoff))
            for i in size_to_aa: 
                try:
                    net = nx.read_gpickle(jn(subfolder, str(cutoff)+'_'+str(thresh)+'.p'))
                    nodes_studied = [node for node in net.nodes if node[0] in size_to_aa[i]]
                    WK[num][cutoff-3, i] = np.average(np.divide(np.asarray(net.degree(nodes_studied, weight='weight')).transpose()[1].astype(np.float),np.asarray(net.degree(nodes_studied)).transpose()[1].astype(np.float)))
                except:
                    WK[num][cutoff-3, i] = None
        
        def plot_avg(M, name):
            f = plt.figure()
            lines = []
            for i in size_to_aa:
                color = next(colors)
                if i%7==0:
                    linestyle=next(style)
                plt.plot(range(3,10), M[:,i], color=color, linestyle=linestyle)
                lines.append(mlines.Line2D([], [], color=color, label=i, linestyle=linestyle))
                f.legend(handles=lines, loc='best', fontsize=8, bbox_to_anchor=(0.5, 1), ncol=4, fancybox=True, shadow=True)
                plt.xlabel('Cutoff distance (in $\AA)$')
            plt.savefig(jn(OUT_FOLDER, name))
        
        plot_avg(WK[num], 'WK_aa_'+folder+'.svg') 
        color=next(colors)
    plot_avg(WK[1]-WK[0], 'WK_aa_diff.svg')
