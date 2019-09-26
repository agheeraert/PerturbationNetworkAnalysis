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
apo_folder = jn(BASE_FOLDER, 'apo')
prfar_folder = jn(BASE_FOLDER, 'prfar')
OUT_FOLDER = '/home/agheerae/results/local_nw/noH/'
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))


def nw(net, node):
    return net.degree(node)/net.degree(node, weight='weight')

def delta_nw(net2, net1, node):
    return nw(net2, node)-nw(net1,node)

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
                new_weight = nw(apo_net, u)+nw(apo_net,v)
            else:
                new_weight = nw(prfar_net, u)+nw(prfar_net,v)
            if new_weight > 0:
                new_net.add_edge(u, v, color='red', weight=abs(new_weight))
            elif new_weight < 0:
                new_net.add_edge(u, v, color='blue', weight=abs(new_weight))
        nx.write_gpickle(new_net, jn(OUT_FOLDER, str(cutoff)+'_0.p'))
        empty, thresh = False, 1
        while not empty:
            L_edge_remove = []
            for u, v in new_net.edges():
                if new_net.get_edge_data(u, v)['weight'] < thresh:
                    L_edge_remove.append((u,v))
            new_net.remove_edges_from(L_edge_remove)
            new_net.remove_nodes_from(list(nx.isolates(new_net)))
            if len(new_net.nodes()) == 0:
                empty = True
            else:
                nx.write_gpickle(new_net, jn(OUT_FOLDER, str(cutoff)+'_'+str(thresh)+'.p'))
            thresh += 1 
