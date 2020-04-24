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
from mayavi import mlab
from tqdm import tqdm
import warnings
import argparse
from copy import deepcopy
from CreateNetwork import three2one
from time import sleep
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

mlab.options.offscreen = True

class DrawNetwork():
    """Draws (and save) perturbation networks"""
    def __init__(self, network, output, pdb_path=None, method='default', colors=['red', 'dodgerblue'], single=False, increment=1, node_size=100, font_size=10):
        print('Drawing...')
        self.net = network
        self.colors = colors
        self.pdb_path = pdb_path
        self.font_size = font_size
        self.node_size = node_size
        self.output = output
        self.single = single
        self.increment = increment
        if method == 'IGPS':
            self.draw_IGPS()
        elif method == '4CFF':
            self.draw_4CFF()
        elif method == '3d':
            assert self.pdb_path, 'No PDB file specified'
            self.draw_default(pos=None, threeD=True)
        elif method == 'single':
            self.drawing(self.output, pos=None)
        else:
            self.draw_default(pos=None)

    def draw_4CFF(self):
        """A specific drawing adapted for 4CFF"""
        structure = PDBParser().get_structure('X', self.pdb_path)[0]
        pos = {}
        for atom in structure.get_atoms():
            if atom.id == 'CA':
                residue = atom.parent
                if residue.resname in three2one:
                    db = sqrt((atom.coord[0]-37.726)**2+(atom.coord[1]-23.815)**2+(atom.coord[2]-34.728)**2)
                    dc = sqrt((atom.coord[0]-32.102)**2+(atom.coord[1]-24.797)**2+(atom.coord[2]-59.026)**2)
                    dd = sqrt((atom.coord[0]-70.621)**2+(atom.coord[1]-72.942)**2+(atom.coord[2]-43.637)**2)
                    y = sqrt(abs(dc**2-((dc**2-db**2)/(2*24.959)+24.959/2)**2))
                    x = sqrt(abs(dd**2-((dd**2-db**2)/(2*59.79)+59.79/2)**2))
                    pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)
        if not self.single:
            self.draw_default(pos)
        else:
            self.drawing(self.output, pos)
       

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
        if not self.single:       
            self.draw_default(pos)
        else:
            self.drawing(self.output, pos)   

    def draw_default(self, pos, threeD=False):
        """Default method to draw the networks. A dictionnary (pos) can be given to draw specific nodes in specific positions"""
        threshold = 0
        empty = False
        while not empty:
            to_remove_edges = []
            for u, v in self.net.edges():
                weight = self.net.get_edge_data(u, v)['weight']
                if weight < threshold:
                    to_remove_edges.append((u, v))
            self.net.remove_edges_from(to_remove_edges)
            self.net.remove_nodes_from(list(nx.isolates(self.net)))

            if len(self.net.edges()) != 0:
                output_path = jn(self.output, str(threshold))
                if not threeD:
                    self.drawing(output_path, pos)
                else:
                    pos = self.get_3d_pos()
                    self.draw_3d(output_path, pos)    
                threshold += self.increment 
            else:
                empty = True
            
    def drawing(self, output_path, pos):
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
            nx.draw(self.net, font_weight='bold', labels={node: node[:-2] for node in self.net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=self.node_size, node_shape='o', font_size=self.font_size, node_color='lightgrey')
        else:
            nx.draw(self.net, font_weight='bold', labels={node: node[:-2] for node in self.net.nodes()}, width=width, edge_color=colors, node_size=self.node_size, node_shape='o', font_size=self.font_size, node_color='lightgrey')
        nx.write_gpickle(self.net, output_path+'.p')
        plt.savefig(output_path+'.pdf')
        plt.close()

    def draw_3d(self, output_path, pos):
        mlab.clf()
        figure = mlab.figure(bgcolor=(1, 1, 1))
        name2id = {elt: i for i, elt in enumerate(self.net.nodes())}
        id2name = {i: elt for i, elt in enumerate(self.net.nodes())}
        integer_net = nx.relabel_nodes(self.net, name2id)
        xyz = np.array([pos[elt] for elt in self.net.nodes()])
        figure.scene.disable_render = True # Super duper trick
        pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                    color=(0.5,0.5,0.5),
                    scale_factor=1,
                    scale_mode='none',
                    resolution=20)
        for i, x in enumerate(xyz):
            mlab.text3d(x[0], x[1], x[2], id2name[i], scale=(1, 1, 1), color=(0,0,0))
        colors_str = list(nx.get_edge_attributes(integer_net, 'color').values())
        color_str2rgb = {'g': (0., 0.5, 1.), 'r': (1., 0., 0.)}
        widths = np.array(list(nx.get_edge_attributes(integer_net, 'weight').values()))
        colors_rgb = [color_str2rgb[color] for color in colors_str]

        # pts.mlab_source.dataset.lines = np.array(integer_net.edges())
        # tube = mlab.pipeline.tube(pts, tube_radius=0.1, tube_color=colors_rgb)
        # mlab.pipeline.surface(tube)
        for i, (u, v) in enumerate(tqdm(integer_net.edges())):
            mlab.plot3d(*[[xyz[u, i], xyz[v, i]] for i in range(0,3)], tube_radius=widths[i]*0.05, color=color_str2rgb[colors_str[i]])
        mlab.savefig(output_path+'.png', magnification=10)
        nx.write_gpickle(self.net, output_path+'.p')

    def get_3d_pos(self):
        structure = PDBParser().get_structure('X', self.pdb_path)[0]
        pos = {}        
        for atom in structure.get_atoms():
            if atom.id == 'CA':
                residue = atom.parent
                if residue.resname in three2one:
                    pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = atom.coord
        return pos

