import os
from os import listdir
from os.path import join as jn
from os.path import dirname, basename
import networkx as nx
import argparse
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import matplotlib
import pickle as pkl

from DrawNetwork import DrawNetwork
from CreateNetwork import three2one
from CreatePerturbationNetwork import CreatePerturbationNetwork

parser = argparse.ArgumentParser(description='Create induced perturbation network from perturbation network.')
parser.add_argument('-f',  type=str, nargs='+',
                    help='List of network to redraw with the zoomed induced perturbation network')
parser.add_argument('-drawing_method',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-pdb_path',  type=str,
                    help='PDB structure file to help draw the network (works only with the IGPS drawing method as of now)')
parser.add_argument('-drawing_colors',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')
args = parser.parse_args()

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

colors = itertools.cycle(('orange', 'g', 'b', 'magenta', 'r', 'black'))   # 'purple',  between olive and navy
tuple_residues = [(['R', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y', 'H'], 'p'),
                (['I', 'L', 'V', 'M', 'F', 'W', 'C', 'P', 'G', 'A'], 'h')]

dict_interactions = {('+', '-'): 'salt_bridge',
                    ('-', '+'): 'salt_bridge',
                    ('+', '+'): 'same_charge',
                    ('-', '-'): 'same_charge',
                    ('h', 'h'): 'hydrophobic',
                    ('p', 'p'): 'polar'
}

dict_residues = {}
for tup in tuple_residues:
    for elt in tup[0]:
        dict_residues[elt] = tup[1]

for fichier in args.f:
    wd = dirname(fichier)
    charged_net = nx.Graph()
    for type_of_interaction in dict_interactions:
        out_folder = jn(wd, dict_interactions[type_of_interaction])
        mkdir(out_folder)
        net = nx.read_gpickle(fichier)
        to_remove_edges = []
        for u, v in net.edges():
            if (dict_residues[u[0]], dict_residues[v[0]]) in dict_interactions:
                toi = dict_interactions[(dict_residues[u[0]], dict_residues[v[0]])]
                if toi == dict_interactions[type_of_interaction] or toi[::-1] == dict_interactions[type_of_interaction]:
                    pass 
                else:
                    to_remove_edges.append((u, v))
            else:
                to_remove_edges.append((u, v))
        net.remove_edges_from(to_remove_edges)
        net.remove_nodes_from(list(nx.isolates(net)))
        if dict_interactions[type_of_interaction] in ['salt_bridge', 'same_charge']:
            charged_net = nx.compose(charged_net, net)
        if len(net.nodes()) != 0:
            DrawNetwork(net, jn(out_folder, basename(fichier).split('.')[0]), pdb_path=args.pdb_path, method=args.drawing_method, colors=args.drawing_colors, single=True)
    mkdir(jn(wd, 'charged_network'))
    DrawNetwork(charged_net, jn(wd, 'charged_network', basename(fichier).split('.')[0]), pdb_path=args.pdb_path, method=args.drawing_method, colors=args.drawing_colors, single=True)
