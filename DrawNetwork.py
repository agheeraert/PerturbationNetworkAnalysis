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

class DrawNetwork():
    """Draws (and save) perturbation networks"""
    def __init__(self, network, output, pdb_path=None, method='default', colors=['red', 'dodgerblue']):
        self.net = network
        self.colors = colors
        self.pdb_path = pdb_path
        self.output = output
        if method == 'default':
            self.draw_default(pos=None)
        if method == 'IGPS':
            self.draw_IGPS()
    
    def draw_IGPS(self):
        """A specific drawing adapted to IGPS this should be tunable for other proteins"""
        structure = PDBParser().get_structure('X', self.pdb_path)[0]
        pos = {}
        distance_thresh = 1
        for atom in structure.get_atoms():
            if atom.id == 'CA':
                residue = atom.parent
                c = 1*(residue.parent.id == 'H') #Separate hisH and hisF
                if residue.resname in three2one:
                        "these values represents the 2D projection of IGPS in our classical view"
                        y = (0.1822020302*atom.coord[0] + 0.6987674421*atom.coord[1] - 0.6917560857*atom.coord[2])*(1-0.5*c)
                        x = 0.9980297273*atom.coord[0]+ 0.0236149631*atom.coord[1]+ 0.05812914*atom.coord[2]
                        pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)
        self.draw_default(pos)
   
    def draw_default(self, pos):
        """Default method to draw the networks. A dictionnary (pos) can be given to draw specific nodes in specific positions"""
        threshold = 0
        empty = False
        print('wesh')
        while not empty:
            to_remove_edges = []
            for u, v in self.net.edges():
                weight = self.net.get_edge_data(u, v)['weight']
                if weight < threshold:
                    to_remove_edges.append((u, v))
            self.net.remove_edges_from(to_remove_edges)
            self.net.remove_nodes_from(list(nx.isolates(self.net)))

            if len(self.net.edges()) != 0:
                print('ca trace')
                fig = plt.figure()
                colors = list(nx.get_edge_attributes(self.net, 'color').values())
                for i, color in enumerate(colors):
                    if color == 'g':
                        colors[i] = self.colors[1]
                    if color == 'r':
                        colors[i] = self.colors[0]                
                width = list(nx.get_edge_attributes(self.net, 'weight').values())
                max_width = max(width)
                for i, elt in enumerate(width):
                    width[i] = elt/max_width*5
                weights = {(u, v): round(nx.get_edge_attributes(self.net, 'weight')[(u,v)]) for (u, v) in nx.get_edge_attributes(self.net, 'weight')}
                if pos:
                    pos = {node: pos[node] for node in self.net.nodes()}
                    nx.draw(self.net, font_weight='bold', nodelist=[node for node in self.net.nodes() if node[-1]=='F'], labels={node: node[:-2] for node in self.net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')
                    nx.draw(self.net, font_weight='bold', nodelist=[node for node in self.net.nodes() if node[-1]=='H'], labels={node: node[:-2] for node in self.net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')
                else:
                    nx.draw(self.net, font_weight='bold', nodelist=[node for node in self.net.nodes() if node[-1]=='F'], labels={node: node[:-2] for node in self.net.nodes()}, width=width, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')
                    nx.draw(self.net, font_weight='bold', nodelist=[node for node in self.net.nodes() if node[-1]=='H'], labels={node: node[:-2] for node in self.net.nodes()}, width=width, edge_color=colors, node_size=100, node_shape='o', font_size=10, node_color='lightgrey')                        
                output_path = jn(self.output, str(threshold)) 
                nx.write_gpickle(self.net, output_path+'.p')
                plt.savefig(output_path+'.pdf')
                threshold+=1
                plt.close()
            else:
                empty = True            

