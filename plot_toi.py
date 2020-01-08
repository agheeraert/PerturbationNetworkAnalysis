from os import listdir
import os
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import matplotlib
import pickle as pkl
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Plot the type of interactions in function of the cutoff')
parser.add_argument('-f',  type=str, nargs='+',
                    help='List of input folders')
parser.add_argument('-o',  type=str,
                    help='Output folder')
parser.add_argument('-s',  type=str, nargs='+',
                    help='List of strings related to each input folder')                    
args = parser.parse_args()

OUT_FOLDER = args.o
def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
mkdir(OUT_FOLDER)

colors = itertools.cycle(('orange', 'g', 'b', 'magenta', 'r', 'black'))   # 'purple',  between olive and navy
tuple_residues = [(['R', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y', 'H'], 'p'),
                (['I', 'L', 'V', 'M', 'F', 'W', 'C', 'P', 'G', 'A'], 'h')]
aa_list = [aa for tup in tuple_residues for aa in tup[0]]
q = len(aa_list)
aa_to_id = dict(zip(aa_list, range(q)))

dict_residues = {}
for tup in tuple_residues:
    for elt in tup[0]:
        dict_residues[elt] = tup[1]

def update_dico(dico, tup):
    if tup in dico:
        dico[tup]+=1
    elif tup[::-1] in dico:
        dico[tup[::-1]]+=1
    else:
        dico[tup] = 1
    return dico

dict_interactions = {('+', '-'): 'sb',
                    ('+', '+'): 'un_sb',
                    ('-', '-'): 'un_sb',
                    ('h', 'h'): 'h',
                    ('p', 'p'): 'p'
}
# matplotlib.rcParams['font.family'] = "serif"
if __name__ =='__main__':
    for i, folder in enumerate(args.f):
        L_threshs = sorted([int(file.split('.')[0]) for file in listdir(folder) if file[-1]=='p'])
        interactions_count = {'sb': ((len(L_threshs))*[0], "Different charge (salt bridge)"),
            'un_sb': ((len(L_threshs))*[0], "Same charge"),
            'h': ((len(L_threshs))*[0], "Hydrophobic"),
            'p': ((len(L_threshs))*[0], "Polar"),
            '?': ((len(L_threshs))*[0], "Non-specific"),
            'c': ((len(L_threshs))*[0], "Covalently bound"),
                    }
        c_tot = (len(L_threshs))*[0]
        for threshold in L_threshs:
            thresh_str = str(threshold)
            if thresh_str+'.p' in listdir(folder):
                net = nx.read_gpickle(jn(folder, thresh_str+'.p'))
                for u, v in net.edges():
                    if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) == 1:
                        interaction = 'c'
                    else:
                        residue_residue = (dict_residues[u[0]], dict_residues[v[0]])
                        if residue_residue in dict_interactions:
                            interaction = dict_interactions[residue_residue]
                        elif residue_residue[::-1] in dict_interactions:
                            interaction = dict_interactions[residue_residue[::-1]]
                        else:
                            interaction = '?'
                    interactions_count[interaction][0][threshold]+=1
                    c_tot[threshold]+=1
            for elt in interactions_count:
                if c_tot[threshold] != 0:
                    interactions_count[elt][0][threshold] = interactions_count[elt][0][threshold]/c_tot[threshold]*100
        lines = []
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.tick_params(labelsize='large')
        ax.set_ylim(0, 101)
        for elt in interactions_count:
            color = next(colors)
            lines.append(mlines.Line2D([], [], color=color, label=interactions_count[elt][1]))
            ax.plot(L_threshs, interactions_count[elt][0], c=color)
        ax2 = ax.twinx()
        ax2.tick_params(labelsize='large')
        ax2.semilogy(L_threshs, c_tot, color='black', linestyle='dotted')
        # ax2.set_ylim(1, 100000)
        lines.append(mlines.Line2D([], [], color='black', label='Total', linestyle=':'))
        ax.grid(linestyle='--')
        #f.legend(handles=lines, loc='upper center', fontsize=8, bbox_to_anchor=(0.5, 1),
        #    ncol=2, fancybox=True, shadow=True)
        plt.tight_layout()
        plt.savefig(jn(OUT_FOLDER, 'toi_'+args.s[i]+'.png'))