import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser

three2one = dict(zip(aa3, aa1))
# path = '/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'
# path = '/home/aria/Stage4A_partie1/python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'
path = '/home/aria/PerturbationNetworkAnalysis/data/frame1000/frame_1000_apo.pdb'
output = '/home/aria/PerturbationNetworkAnalysis/data/frame1000/frame_1000.tcl'
A = nx.read_gpickle('tests/test.p')
structure = PDBParser().get_structure('X', path)[0]

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords
with open(output, 'w') as output_file:
        output_file.write('draw color 1 \n')
        for u, v in A.edges():
            if u[-1] != v[-1]:
                output_file.write('draw color 3 \n')
            try:
                output_file.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius 1 \n')
            except KeyError:
                pass
            if u[-1] != v[-1]:
                output_file.write('draw color 1 \n')
                
        for u in A.nodes():
            output_file.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius 1.5 \n')