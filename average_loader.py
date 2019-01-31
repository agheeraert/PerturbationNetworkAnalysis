import networkx as nx
import numpy as np
from os.path import dirname
from os.path import join as jn
import matplotlib.pyplot as plt
import itertools
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt, copysign
import warnings
import argparse
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))
parser = argparse.ArgumentParser(description='Makes the difference between two graph to get their perturbation network')
parser.add_argument('path1',  type=str,
                     help='first file to load')
parser.add_argument('path2',  type=str,
                     help='second file to load')
parser.add_argument('output',  type=str,
                     help='output folder + base string')                               
parser.add_argument('-pdb',  type=str, default='/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/1frame_1.pdb',
                     help='file structure to draw on')                   
args = parser.parse_args()

net1, net2 = nx.read_gpickle(args.path1), nx.read_gpickle(args.path2)
perturbation = nx.compose(net1, net2)
empty = False
threshold = 0
while not empty:
    mapping_weight, mapping_color = {}, {}
    for u, v in set(net2.edges) & set(net1.edges) & set(perturbation.edges):
        w = net2.get_edge_data(u, v)['weight'] - net1.get_edge_data(u,v)['weight']
        if abs(w) > threshold:
            mapping_weight[(u, v)] = abs(w)
            if w > threshold:
                mapping_color[(u, v)] = 'r'
            else:
                mapping_color[(u, v)] = 'g'
        else:
            perturbation.remove_edge(u,v)
    for u, v in (set(net2.edges) - set(net1.edges)) & set(perturbation.edges):
        weight_edge = net2.get_edge_data(u, v)['weight']
        if weight_edge > threshold:
            mapping_weight[(u, v)] = weight_edge
            mapping_color[(u, v)] = 'r'
        else:
            perturbation.remove_edge(u,v)
    for u, v in (set(net1.edges) - set(net2.edges)) & set(perturbation.edges):
        weight_edge = net1.get_edge_data(u, v)['weight']
        if weight_edge > threshold:
            mapping_weight[(u, v)] = weight_edge
            mapping_color[(u, v)] = 'g'
        else:
            perturbation.remove_edge(u,v)
    nx.set_edge_attributes(perturbation, name='color', values=mapping_color)
    nx.set_edge_attributes(perturbation, name='weight', values=mapping_weight)
    perturbation.remove_nodes_from(list(nx.isolates(perturbation)))
    if len(perturbation.nodes()) == 0:
        empty = True
    else:
        nx.write_gpickle(perturbation, args.output+str(threshold)+'.p')
    threshold+=1