import os 
from os.path import dirname, basename
from os.path import join as jn
import argparse
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import networkx as nx

from DrawNetwork import DrawNetwork
from CreateNetwork import three2one
from CreatePerturbationNetwork import PerturbationNetwork

parser = argparse.ArgumentParser(description='Create induced perturbation network from perturbation network.')
parser.add_argument('-f',  type=str, nargs='+',
                    help='List of network to redraw with the zoomed induced perturbation network')
parser.add_argument('-r', type=str, nargs='+',
                    help='Root nodes used for the induced perturbation network')
parser.add_argument('-dm',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-p',  type=str,
                    help='PDB structure file to help draw the network')
parser.add_argument('-dc',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')
parser.add_argument('-fs',  type=int, default=10,
                    help='Font size')
parser.add_argument('-ns',  type=int, default=100,
                    help='Node size')
args = parser.parse_args()

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for _file in args.f:
    net = nx.read_gpickle(_file)
    _weights = nx.get_edge_attributes(net, 'weight')
    _colors = nx.get_edge_attributes(net, 'color')
    for root in args.r:
        if root in net.nodes():
            out_folder = jn(dirname(_file), root)
            mkdir(out_folder)
            tree = nx.Graph(nx.bfs_tree(net, root))
            for u, v in net.edges():
                if u in tree.nodes() and v in tree.nodes() and not (u, v) in tree.edges() and not (v, u) in tree.edges():
                    tree.add_edge(u, v)
            nx.set_edge_attributes(tree, name='weight', values=_weights)
            nx.set_edge_attributes(tree, name='color', values=_colors)
            DrawNetwork(tree, jn(out_folder, basename(_file).split('.')[0]), pdb_path=args.p, method=args.dm, colors=args.dc, single=True, node_size=args.ns, font_size=args.fs)
