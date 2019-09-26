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

BASE_FOLDER = '/home/agheerae/results/avg_sims/pertnet/'
OUT_FOLDER = '/home/agheerae/results/local_nw/noH/'
endword = 'H'

q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))

def nw(node):
    return node.degree()/node.degree(weight='weight')

if __name__ =='__main__':
    W, K, WK = np.zeros((7, 6)), np.zeros((7, 6)), np.zeros((7, 6))
    for cutoff in range(3, 10):
        folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
        for thresh in range(1, 12, 2): 
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+str(thresh)+'.p'))
            new_net = nx.Graph()
            for u, v in net.edges():
                new_weight=abs(nw(u)+nw(v))
                if nw(u)+nw(u)>0:
                    new_net.add_edge(u, v, color='red', weight=new_weight)
                elif nw(u)+nw(v)<0:
                    new_net.add_edge(u, v, color='blue', weight=new_weight)
            nx.write_gpickle(new_net, jn(OUT_FOLDER, str(cutoff)+'_'+str(thresh)+'.p'))



