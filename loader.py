import networkx as nx
import numpy as np
import os
from os.path import dirname
from os.path import join as jn
import matplotlib.pyplot as plt
import itertools
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt, copysign
import warnings
import argparse
from copy import deepcopy
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))

parser = argparse.ArgumentParser(description='Load graphs and rebuild them')
parser.add_argument('path',  type=str,
                     help='file to load')
parser.add_argument('-pdb',  type=str, default='/home/agheerae/Article/base_IGPS.pdb',
                     help='file structure to draw on')
parser.add_argument('-r',  type=bool, default=False,
                     help='Indicates there\'s a root in the filename')                    
parser.add_argument('-t',  type=int , default=None,
                     help='Splits the 2D graph representation until a specific certain threshold')
args = parser.parse_args()

def split_man(pos, L_nodes, L_direction_x, L_direction_y):
    for i, node in enumerate(L_nodes):
        pos[node] = (pos[node][0]+L_direction_x[i], pos[node][1]+L_direction_y[i])
    return pos

def split(pos, d_thresh, n_max=100):
    n = 0
    conflicts = ['a', 'b']
    while len(conflicts) >= 1 and n < n_max: 
        conflicts = []
        L_to_move = []
        for elt1 in pos:
            for elt2 in pos:
                if elt1 != elt2:
                    if elt1 not in L_to_move and elt2 not in L_to_move:
                        if abs(pos[elt1][1] - pos[elt2][1]) < d_thresh and abs(pos[elt1][0] - pos[elt2][0]) < d_thresh:
                            conflicts.append((elt1, elt2))
                            L_to_move.extend([elt1, elt2])
        for elt1, elt2 in conflicts:
            a = ((pos[elt1][0] - pos[elt2][0]) >= 0)*1 
            b = ((pos[elt1][1] - pos[elt2][1]) >= 0)*1 
            pos[elt1] = (pos[elt1][0] + a*0.1, pos[elt1][1] + b*0.1)
            pos[elt2] = (pos[elt2][0] - a*0.1, pos[elt2][1] - b*0.1)
        n+=1
    return pos

def split2(pos, d_thresh, nmax=500):
    n = 0
    n_conf = 1
    while n < nmax and n_conf != 0:
        n_conf = 0
        for elt1 in pos:
            for elt2 in pos:
                if elt1 != elt2:
                    if sqrt((pos[elt1][1]-pos[elt2][1])**2+(pos[elt1][0]-pos[elt2][0])**2) < d_thresh:
                        dx_sign = (pos[elt1][0]-pos[elt2][0] > 0)*1
                        pos[elt1] = (pos[elt1][0] + dx_sign*d_thresh/2, pos[elt1][1])
                        pos[elt2] = (pos[elt2][0] - dx_sign*d_thresh/2, pos[elt2][1])
                        n_conf+=1
                    if sqrt((pos[elt1][1]-pos[elt2][1])**2+(pos[elt1][0]-pos[elt2][0])**2) < d_thresh:
                        dy_sign = (pos[elt1][1]-pos[elt2][1] > 0)*1
                        pos[elt1] = (pos[elt1][0], pos[elt1][1]+dy_sign*d_thresh/2)
                        pos[elt2] = (pos[elt2][0], pos[elt2][1]-dy_sign*d_thresh/2)
    return pos

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

structure = PDBParser().get_structure('X', args.pdb)[0]
pos = {}
distance_thresh = 1
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        c = 1*(residue.parent.id == 'H')
        print(residue.id[1], residue.resname)
        if residue.resname in three2one:
                y = (0.1822020302*atom.coord[0] + 0.6987674421*atom.coord[1] - 0.6917560857*atom.coord[2])*(1-0.5*c)
                x = 0.9980297273*atom.coord[0]+ 0.0236149631*atom.coord[1]+ 0.05812914*atom.coord[2]
                pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)

output_str1 = args.path.split('_')[0]
if args.r:
    root = args.path.split('_')[2].split('.')[0]
net = nx.read_gpickle(args.path)
output_folder = dirname(args.path)
threshold = 1
empty = False
while not empty:
    to_remove_edges = []
    for u, v in net.edges():
        weight = net.get_edge_data(u, v)['weight']
        if weight < threshold:
            to_remove_edges.append((u, v))
    net.remove_edges_from(to_remove_edges)
    net.remove_nodes_from(list(nx.isolates(net)))

    if len(net.edges()) != 0:
        fig = plt.figure(figsize=[6.4, 4.8], dpi=200)
        colors = list(nx.get_edge_attributes(net, 'color').values())
        for i, color in enumerate(colors):
            if color == 'g':
                colors[i] = 'dodgerblue'                
        width = list(nx.get_edge_attributes(net, 'weight').values())
        max_width = max(width)
        for i, elt in enumerate(width):
            width[i] = elt/max_width*5
        weights = {(u, v): round(nx.get_edge_attributes(net, 'weight')[(u,v)]) for (u, v) in nx.get_edge_attributes(net, 'weight')}
        pos = {node: pos[node] for node in net.nodes()}
        if threshold == args.t:
            pos = split2(pos, 2.5)
            pos['R18:H'] = pos['Y17:H']
            pos = split_man(pos, ['Y182:F', 'N103:F', 'V12:F', 'V18:F', 'F227:F', 'Y17:H', 'L42:H', 'L118:H'], [-2, -2, 4, 2, 0, 0, 0, 0], [0, 0, 4, 0, 2, 2, -1, 1])
        nx.draw(net, font_weight='bold', nodelist=[node for node in net.nodes() if node[-1]=='F'], labels={node: node[:-2] for node in net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')
        # nx.draw_networkx_edge_labels(net, pos=pos, edge_labels=weights, font_color='black', font_size=5)
        nx.draw(net, font_weight='bold', nodelist=[node for node in net.nodes() if node[-1]=='H'], labels={node: node[:-2] for node in net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')
        thresh = str(round(threshold, 1)).replace('.', '-')
        if not args.r:
            output = output_str1+'_'+thresh
        else:
            output = output_str1+'_'+thresh+'_'+root
        nx.write_gpickle(net, output+'.p')
        plt.savefig(output+'.pdf')
        threshold+=1
    else:
        empty = True



