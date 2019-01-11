import networkx as nx
import numpy as np
from os.path import dirname
import matplotlib.pyplot as plt
import itertools
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))

pdb = '/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/1frame_1.pdb'
structure = PDBParser().get_structure('X', pdb)[0]
pos = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                y = atom.coord[1]
                x = sqrt(2)*atom.coord[2] - sqrt(2)*atom.coord[0]
                pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)

for i in range(9, 10):
    load_path = '/home/agheerae/results/sim4/pertnet/cutoff_'+str(i)+'/'+str(i)+'_0.p'
    output_str1 = '/home/agheerae/results/sim4/pertnet/cutoff_'+str(i)+'/'+str(i)+'_'
    net = nx.read_gpickle(load_path)
    output_folder = dirname(load_path)
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
            weights = nx.get_edge_attributes(net, 'weight')
            for weight in weights:
                weights[weight] = str(weights[weight])
            nx.draw(net, with_labels=True, font_weight='bold', width=width, pos=pos, edge_labels=weights, edge_color=colors, node_size=100, node_color='grey', font_size=8)
            thresh = str(round(threshold, 1)).replace('.', '-')
            nx.write_gpickle(net, output_str1+thresh+'.p')
            plt.savefig(output_str1+thresh+'.svg')

            threshold+=1
        else:
            empty = True



