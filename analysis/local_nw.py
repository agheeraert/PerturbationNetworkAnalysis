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

BASE_FOLDER = '/home/agheerae/results/avg_sims/nets/'
APO_FOLDER = jn(BASE_FOLDER, 'apo')
PRFAR_FOLDER = jn(BASE_FOLDER, 'prfar')
OUT_FOLDER = '/home/agheerae/results/local_nw/noH/'
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))

def delta_nw(net2, net1, node):
    return net2.degree(node)/net2.degree(node, weight='weight')-net1.degree(node)/net1.degree(node, weight='weight')

if __name__ =='__main__':
    W, K, WK = np.zeros((7, 6)), np.zeros((7, 6)), np.zeros((7, 6))
    for cutoff in range(3, 10):
        apo_net = nx.read_gpickle(jn(apo_folder, str(cutoff)+'.p'))
        prfar_net = nx.read_gpickle(jn(prfar_folder, str(cutoff)+'.p'))
        new_net = nx.Graph()
        for u, v in nx.compose(apo_net, prfar_net).edges():
            if (u,v) in apo_net and (u, v) in prfar_net:
                new_weight=delta_nw(prfar_net, apo_net, u)+nw(prfar_net, apo_net, net, v)
            elif (u,v) in apo_net and (u,v) not in prfar_net:
                new_weight = delta_nw(apo_net)
            else:
                new_weight = delta_nw(prfar_net)
            if new_weight > 0:
                new_net.add_edge(u, v, color='red', weight=abs(new_weight))
            elif new_weight < 0:
                new_net.add_edge(u, v, color='blue', weight=abs(new_weight))
        nx.write_gpickle(new_net, jn(OUT_FOLDER, str(cutoff)+'.p'))



