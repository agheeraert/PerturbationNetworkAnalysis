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
from Bio.PDB.Polypeptide import aa3
import argparse

parser = argparse.ArgumentParser(description='Plot the average weight, degree and nw variation of the cutoff')
parser.add_argument('f',  type=str,
                    help='Input folder')
parser.add_argument('o',  type=str,
                    help='Output path')  
parser.add_argument('-s',  type=str, default='',
                    help='Output string')                 
args = parser.parse_args()


colors = itertools.cycle(('orange', 'g', 'b', 'magenta', 'r', 'black'))   # 'purple',  between olive and navy
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))

def plot_avg(M, name):
    f = plt.figure()
    lines = []
    for thresh in range(1,12,2):
        color = next(colors)
        plt.plot(range(3,10), M.transpose()[thresh//2], color=color)
        lines.append(mlines.Line2D([], [], color=color, label='$\\bar{w}$='+str(thresh)))
        f.legend(handles=lines, loc='upper center', fontsize=8, bbox_to_anchor=(0.5, 1), ncol=2, fancybox=True, shadow=True)
        plt.xlabel('Cutoff distance (in $\AA)$')
        plt.ylabel('Average '+name.split('.')[0].lower()+' per node in the perturbation network')
    plt.savefig(jn(args.o, name))


if __name__ =='__main__':
    W, K, WK = np.zeros((7, 6)), np.zeros((7, 6)), np.zeros((7, 6))
    for cutoff in range(3, 10):
        folder = jn(args.f, 'cutoff'+str(cutoff))
        for thresh in range(1, 12, 2): 
            try:
                net = nx.read_gpickle(jn(folder, str(thresh)+'.p'))
                K[cutoff-3,(thresh//2)] = np.average(np.asarray(net.degree(net.nodes)).transpose()[1].astype(np.float))
#                W[cutoff-3,(thresh//2)] = np.average(np.asarray([net.get_edge_data(u, v)['weight'] for u, v in net.edges()]))
                W[cutoff-3,(thresh//2)] = np.average(np.asarray(net.degree(net.nodes, weight='weight')).transpose()[1].astype(np.float))
                WK[cutoff-3,(thresh//2)] = np.average(np.divide(np.asarray(net.degree(net.nodes, weight='weight')).transpose()[1].astype(np.float),np.asarray(net.degree(net.nodes)).transpose()[1].astype(np.float)))
            except:
                K[cutoff-3,(thresh//2)] = None
                W[cutoff-3,(thresh//2)] = None
                WK[cutoff-3,(thresh//2)] = None   
    plot_avg(K, 'K'+args.s+'.png') 
    plot_avg(W, 'W'+args.s+'.png') 
    plot_avg(WK, 'WK'+args.s+'.png') 

