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
OUT_FOLDER = '/home/agheerae/results/individual_mesures_vs_cutoff'
q = len(aa3)
aa_to_id = dict(zip(aa3, range(q)))


def nw(net, node):
    return net.degree(node)/net.degree(node, weight='weight')

def dictup(dic, node, value):
    if node in dic:
        dic[node].append(value)
    else:
        dic[node] = [value]

def plotdic(dico, string):
    plt.figure()
    for i, elt in enumerate(dico):
        plt.subplot(46,10, i+1)
        plt.plot(cutoffs, dico[elt])
        plt.title(elt)
    plt.savefig(jn(OUT_FOLDER, string))
    plt.close()

def deltadic(dic1, dic2):
    dic = {}
    for elt in dic2:
        dic[elt] = (np.asarray(dic2[elt])-np.asarray(dic1[elt])).tolist()
    return dic

if __name__ =='__main__':
    degree_apo = {}
    degree_prfar = {}
    weight_apo = {}
    weight_prfar = {}
    nw_apo = {}
    nw_prfar = {}
    cutoffs = range(3, 10)
    for cutoff in cutoffs:
            apo_net = nx.read_gpickle(jn(apo_folder, str(cutoff)+'.p'))
            prfar_net = nx.read_gpickle(jn(prfar_folder, str(cutoff)+'.p'))
            for node in apo_net.nodes():
                dictup(degree_apo, node, apo_net.degree(node))
                dictup(weight_apo, node, apo_net.degree(node, weight='weight'))
                dictup(nw_apo, node, nw(apo_net, node))
            for node in prfar_net.nodes():
                dictup(degree_prfar, node, prfar_net.degree(node))
                dictup(weight_prfar, node, prfar_net.degree(node, weight='weight'))
                dictup(nw_prfar, node, nw(prfar_net, node))
    degree_delta = deltadic(degree_apo, degree_prfar)
    weight_delta = deltadic(weight_apo, weight_prfar)
    nw_delta = deltadic(nw_apo, nw_prfar)


    plotdic(degree_apo, 'degree_apo.png')
    plotdic(degree_prfar, 'degree_prfar.png')
    plotdic(degree_delta, 'degree_delta.png')
    plotdic(weight_apo, 'weight_apo.png')
    plotdic(weight_prfar, 'weight_prfar.png')
    plotdic(weight_delta, 'weight_delta.png')    
    plotdic(nw_apo, 'nw_apo.png')
    plotdic(nw_prfar, 'nw_prfar.png')
    plotdic(nw_delta, 'nw_delta.png')
        



