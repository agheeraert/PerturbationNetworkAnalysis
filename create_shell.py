import os 
from os.path import join as jn
import networkx as nx
import argparse
import numpy as np
from DrawNetwork import DrawNetwork

parser = argparse.ArgumentParser(description='Create the shell perturbation network difference from two whole perturbation networks')
parser.add_argument('path1',  type=str,
                    help='Input file with the inside of the shell')
parser.add_argument('path2',  type=str,
                    help='Input file with the whole sphere')
parser.add_argument('output',  type=str,
                    help='Folder where to put the shell result')  
parser.add_argument('-drawing_method',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-pdb_path',  type=str,
                    help='PDB structure file to help draw the network')
parser.add_argument('-drawing_colors',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')


args = parser.parse_args()

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

mkdir(args.output)
net1 = nx.read_gpickle(args.path1)
net2 = nx.read_gpickle(args.path2)
outnet = nx.Graph()

def get_relative_weight(u, v, net):
    return (int(net.get_edge_data(u, v)['color']=='r')*2-1)*net.get_edge_data(u, v)['weight']


for u,v in net2.edges():
    if (u, v) in net1.edges():
        new_weight = get_relative_weight(u, v, net2) - get_relative_weight(u, v, net1)
    else:
        new_weight = get_relative_weight(u, v, net2)
    if new_weight > 0:
        new_color = 'r'
    else:
        new_color = 'g'
    if new_weight !=0:
        outnet.add_edge(u, v, weight=abs(new_weight), color=(new_color))


DrawNetwork(outnet, args.output, pdb_path=args.pdb_path, method=args.drawing_method, colors=args.drawing_colors)






