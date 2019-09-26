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
floaters = True
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))

def nw(net, node):
    return net.degree(node)/net.degree(node, weight='weight')

if __name__ =='__main__':
    W, K, WK = np.zeros((7, 6)), np.zeros((7, 6)), np.zeros((7, 6))
    for cutoff in range(3, 10):
        folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
        if cutoff == 3 and floaters: #FLAG3
            L_threshs = sorted([float(file.split('_')[1].split('.')[0].replace('-','.')) for file in listdir(folder) if file[-1]=='p'])
        else:
            L_threshs = sorted([int(file.split('_')[1].split('.')[0]) for file in listdir(folder) if file[-1]=='p'])
 
        for thresh in L_threshs: 
            if cutoff == 3: #FLAG3
                thresh_str = str(thresh).replace('.', '-')
                if thresh == 1.0:
                    thresh_str = '1'
                elif thresh == 0.0:
                    thresh_str = '0'  
                thresh = int(thresh*10)
            else:
                thresh_str = str(thresh)
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+thresh_str+'.p'))
            new_net = nx.Graph()
            for u, v in net.edges():
                new_weight=abs(nw(net, u)+nw(net, v))
                if nw(net, u)+nw(net, u)>0:
                    new_net.add_edge(u, v, color='red', weight=new_weight)
                elif nw(net, u)+nw(net, v)<0:
                    new_net.add_edge(u, v, color='blue', weight=new_weight)
            nx.write_gpickle(new_net, jn(OUT_FOLDER, str(cutoff)+'_'+thresh_str+'.p'))



