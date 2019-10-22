from CreateNetwork import CreateNetwork
from DrawNetwork import DrawNetwork
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
from os.path import join, dirname
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

class CreatePerturbationNetwork(CreateNetwork):
    def __init__(self, path1, path2, pos1=None, pos2=None, cutoff=5, avg=False):
        super().__init__(pos1=pos1, pos2=pos2, cutoff=cutoff)
        if not avg:
            self.net1 = self.create(path1)
            self.net2 = self.create(path2)
        else:
            self.net1 = self.create_average(path1)
            self.net2 = self.create_average(path2)            

    def perturbation(self, threshold):
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
        return perturbation

    def draw_perturbation(self, output, pdb_path=None, method='default', colors=['red', 'dodgerblue']):
        DrawNetwork(self.perturbation(0), output, pdb_path, method=method, colors=['red', 'dodgerblue'])


if __name__ == '__main__':
    CreatePerturbationNetwork(path1='/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/', path2='/home/agheerae/Python/PerturbationNetworkAnalysis/data/prfar_all/', avg=True).draw_perturbation(threshold=range(1, 30, 1), output='/home/agheerae/Python/PerturbationNetworkAnalysis/data/', rearrange=('H', 'F'), save='avg')
