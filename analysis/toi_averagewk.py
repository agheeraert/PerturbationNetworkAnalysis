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
OUT_FOLDER = '/home/agheerae/results/cutoff_analysis/w_and_k/by_amino_acid/'
thresh = 5
colors = itertools.cycle(('orange', 'g', 'b'))   # 'purple',  between olive and navy
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))
dict_residues = { 'c': ['R', 'K', 'D', 'E'],
                  'p': ['S', 'T', 'N', 'Q', 'Y', 'H'],
                  'h': ['I', 'L', 'V', 'M', 'F', 'W', 'C', 'P', 'G', 'A'] }
residues_str = { 'c': 'Charged resides',
                 'p': 'Polar residues',
                 'h': 'Hydrophobic residues'
               }
if __name__ =='__main__':
    WK = [] 
    for num, folder in enumerate(['apo', 'prfar']):
        WK.append(np.zeros((7, 3)))
        for cutoff in range(3, 10):
            subfolder = jn(BASE_FOLDER, folder, 'cutoff_'+str(cutoff))
            for i, typ in enumerate(['c','p','h']): 
                try:
                    net = nx.read_gpickle(jn(subfolder, str(cutoff)+'_'+str(thresh)+'.p'))
                    nodes_studied = [node for node in net.nodes if node[0] in dict_residues[typ]]
                    WK[num][cutoff-3, i] = np.average(np.divide(np.asarray(net.degree(nodes_studied, weight='weight')).transpose()[1].astype(np.float),np.asarray(net.degree(nodes_studied)).transpose()[1].astype(np.float)))
                except:
                    WK[num][cutoff-3, i] = None
        
        def plot_avg(M, name):
            f = plt.figure()
            lines = []
            for i, typ in enumerate(['c','p','h']):
                color = next(colors)
                plt.plot(range(3,10), M[:,i], color=color)
                lines.append(mlines.Line2D([], [], color=color, label=residues_str[typ]))
                f.legend(handles=lines, loc='best', fontsize=8, bbox_to_anchor=(0.5, 1), fancybox=True, shadow=True)
                plt.xlabel('Cutoff distance (in $\AA)$')
            plt.savefig(jn(OUT_FOLDER, name))
        
        plot_avg(WK[num], 'WK_typ_'+folder+'.svg') 
    plot_avg(WK[1]-WK[0], 'WK_typ_diff.svg')

