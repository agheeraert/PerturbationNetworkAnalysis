from os import listdir
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import matplotlib
import pickle as pkl
import pandas as pd
from Bio.PDB.Polypeptide import aa3, aa1
from Bio.PDB import PDBParser
import warnings
from tqdm import tqdm
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import biographs as bg


class CreateLRN():
    def __init__(self, path1, path2, cutoff=5, lrn_threshold=7):
        self.cutoff=cutoff
        self.lrn_threshold = lrn_threshold
        self.three2one = dict(zip(aa3, aa1))
        self.one2three = dict(zip(aa1, aa3))
        self.three2one['5CS'] = 'C'
        self.three2one['HIP'] = 'H'
        self.three2one['HID'] = 'H'
        self.three2one['HIE'] = 'H'
        self.three2one['GLH'] = 'E'
        self.three2one['ASH'] = 'D'
        self.net1 = self.create(path1)
        self.net2 = self.create(path2)

    def remove_close_edges(self, net):
        """ removes edges that are too close in the network """
        L_edge_remove = []
        for u, v in net.edges():
            if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) <= self.lrn_threshold or u[-1] != v[-1]:
                L_edge_remove.append((u,v))
        net.remove_edges_from(L_edge_remove)
        net.remove_nodes_from(list(nx.isolates(net)))
        return net

    def create(self, frame_path):
        """ creates the Long Range Network of a frame """
        mol = bg.Pmolecule(frame_path)
        self.net = mol.network(cutoff=self.cutoff, weight=True)
        self.structure = PDBParser().get_structure('X', frame_path)[0]
        residues = []
        for residue in self.structure.get_residues():
            residues.append(self.three2one[residue.resname])
        old_labels = self.net.nodes
        labels = [a+b[1:]+':'+b[0] for a,b in zip(residues, old_labels)]
        mapping = dict(zip(old_labels, labels))
        self.net = nx.relabel_nodes(self.net, mapping)
        self.net = self.remove_close_edges(self.net)
        return self.net

    def nw(self, net):
        """ returns the nw for all the nodes in a net """
        return {elt[0][0]: elt[0][1]/elt[1][1] for elt in zip(net.degree(weight='weight'), net.degree())}

    def get_nw(self):
        """ get the nw for all the nodes in the built nets"""
        return self.nw(self.net1), self.nw(self.net2)

if __name__ == "__main__":
    for i in range(2,5):
        apo_folder = "/home/aghee/PDB/Apo_frames/Sim"+str(i)
        prfar_folder= "/home/aghee/PDB/Prfar_frames/Sim"+str(i)
        NW1, NW2 = [], []
        for frame1, frame2 in tqdm(zip(sorted(listdir(apo_folder)), sorted(listdir(prfar_folder))), total=1000):
            NWs = CreateLRN(jn(apo_folder, frame1), jn(prfar_folder, frame2)).get_nw()
            NW1.append(NWs[0])
            NW2.append(NWs[1])
        pkl.dump(NW1, open("/home/aghee/results/long_range_evolution/NW_apo_"+str(i)+".p", 'wb'))
        pkl.dump(NW2, open("/home/aghee/results/long_range_evolution/NW_prfar_"+str(i)+".p", 'wb'))



