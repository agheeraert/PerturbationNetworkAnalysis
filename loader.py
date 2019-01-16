import networkx as nx
import numpy as np
from os.path import dirname
from os.path import join as jn
import matplotlib.pyplot as plt
import itertools
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt
import warnings
import argparse
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))

parser = argparse.ArgumentParser(description='Load graphs and rebuild them')
parser.add_argument('path',  type=str,
                     help='file to load')
parser.add_argument('-pdb',  type=str, default='/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/1frame_1.pdb',
                     help='file structure to draw on')
parser.add_argument('-r',  type=bool, default=False,
                     help='Indicates there\'s a root in the filename')                    
args = parser.parse_args()


structure = PDBParser().get_structure('X', args.pdb)[0]
pos = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                y = atom.coord[1]
                x = sqrt(2)*atom.coord[2] - sqrt(2)*atom.coord[0]
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
        fig = plt.figure()
        colors = list(nx.get_edge_attributes(net, 'color').values())
        width = list(nx.get_edge_attributes(net, 'weight').values())
        max_width = max(width)
        for i, elt in enumerate(width):
            width[i] = elt/max_width*2
        weights = {(u, v): round(nx.get_edge_attributes(net, 'weight')[(u,v)], 2) for (u, v) in nx.get_edge_attributes(net, 'weight')}
        nx.draw(net, with_labels=True, font_weight='bold', width=width, pos=pos, edge_color=colors, node_size=100, node_color='grey', font_size=8)
        nx.draw_networkx_edge_labels(net, pos=pos, edge_labels=weights, font_color='black', font_size=5)
        thresh = str(round(threshold, 1)).replace('.', '-')
        if not args.r:
            output = output_str1+'_'+thresh
        else:
            output = output_str1+'_'+thresh+'_'+root
        nx.write_gpickle(net, output+'.p')
        plt.savefig(output+'.svg')
        threshold+=1
    else:
        empty = True



