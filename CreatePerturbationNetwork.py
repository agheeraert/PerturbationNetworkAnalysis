from CreateNetwork import AANetwork
from DrawNetwork import DrawNetwork
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import itertools
import os
from os.path import join, dirname, isdir, basename
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

num = itertools.cycle((1, 2))

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

class PerturbationNetwork(AANetwork):
    """Creates the Perturbation Network"""
    def __init__(self, path1, path2, out_dir, cutoff=5):
        super().__init__(cutoff=cutoff)
        # if not avg:
        #     self.net1 = self.create(path1)
        #     self.net2 = self.create(path2)
        # # elif std:
        # #     self.net1 = self.create_avg_std(path1)
        # #     self.net2 = self.create_avg_std(path2)
        # else:
        #     self.net1 = self.create_average(path1)
        #     self.net2 = self.create_average(path2)
        self.out_dir = out_dir
        def smart_loader(path):
            if len(basename(path)) != 0: #should get the name of a file or a folder without /
                string = basename(path).split('.')[0]+'.p'
            elif len(basename(path[:-1])) !=0: #should get the name of a folder with /
                string = basename(path[:-1])+'.p'
            else: #dummy name generation
                string = 'aa_net%s.p' %next(num)
            if isdir(path): 
                #loads a folder and do the average aanet, saves
                w_dir = join(out_dir, 'aa_net')
                mkdir(w_dir)
                net = self.create_average(path)
                self.save(join(w_dir, string))
            elif path[-4:] == '.pdb': 
                #loads a single frame, do the aanet, saves
                net = self.create(path)
                self.save(join(w_dir, string))
            elif path[-2:] == '.p': 
                #directly loads a precomputed network
                net = nx.read_gpickle(path)
            else:
                print('Unknown extension for file %s' %path)
                return None
            return net
        self.net1 = smart_loader(path1)
        self.net2 = smart_loader(path2)            

    def perturbation(self, threshold=0):
        """Crates the perturbation network between the initialized networks at a specific threshold"""
        perturbation = nx.compose(self.net1, self.net2)
        mapping_weight, mapping_color = {}, {}
        for u, v in set(self.net2.edges) & set(self.net1.edges):
            w = self.net2.get_edge_data(u, v)['weight'] - self.net1.get_edge_data(u,v)['weight']
            if abs(w) > threshold:
                mapping_weight[(u, v)] = abs(w)
                if w > threshold:
                    mapping_color[(u, v)] = 'r'
                else:
                    mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(self.net2.edges) - set(self.net1.edges):
            weight_edge = self.net2.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'r'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(self.net1.edges) - set(self.net2.edges):
            weight_edge = self.net1.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        nx.set_edge_attributes(perturbation, name='color', values=mapping_color)
        nx.set_edge_attributes(perturbation, name='weight', values=mapping_weight)
        perturbation.remove_nodes_from(list(nx.isolates(perturbation)))
        self.perturbation = perturbation

    def draw_perturbation(self, pdb_path=None, method='default', colors=['red', 'dodgerblue']):
        """Draws and save the perturbation network for the initialized networks for all integer 
        thresholds until the network is empty"""
        self.perturbation()
        DrawNetwork(self.perturbation, self.out_dir, pdb_path, method=method, colors=['red', 'dodgerblue'])