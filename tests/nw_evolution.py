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
from Bio.PDB.Polypeptide import aa3
import biographs as bg


class CreateLRN():
    def __init__(self, path1, path2, cutoff=5, lrn_threshold=7):
        self.cutoff=cutoff
        self.lrn_threshold = lrn_threshold
        self.net1 = self.create(path1)
        self.net2 = self.create(path2)

    def remove_close_edges(net):
        """ removes edges that are too close in the network """
        for u, v in net.edges():
            if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) <= self.lrn_threshold or u[-1] != v[-1]:
                L_edge_remove.append((u,v))
        net.remove_edges_from(L_edge_remove)
        net.remove_nodes_from(list(nx.isolates(net)))

    def create(self, frame_path):
        """ creates the Long Range Network of a frame """
        mol = bg.Pmolecule(frame_path)
        self.net = remove_close_edges(mol.network(cutoff=self.cutoff, weight=True))
        self.structure = PDBParser().get_structure('X', pdb)[0]
        residues = []
        for residue in self.structure.get_residues():
            residues.append(self.three2one[residue.resname])
        old_labels = self.net.nodes
        labels = [a+b[1:]+':'+b[0] for a,b in zip(residues, old_labels)]
        mapping = dict(zip(old_labels, labels))
        self.net = nx.relabel_nodes(self.net, mapping)
        return self.net

    def nw(self, net):
        """ returns the nw for all the nodes in a net """
        return net.degree()/net.degree(weight='weight')

    def get_nw(self):
        """ get the nw for all the nodes in the built nets"""
        return nw(self.net1), nw(self.net2)

if __name__ == "__main__":
    apo_folder = "/PDB/Apo_frames/Sim1"
    prfar_folder= "/PDB/Prfar_frames/Sim1"
    NW1 = []
    NW2 = []
    for frame1, frame2 in zip(os.listdir(apo_folder), os.listdir(prfar_folder)):
        NWs = CreateLRN(frame1, frame2).get_nw()
        NW1.append(NW[0])
        NW2.append(NW[1])



