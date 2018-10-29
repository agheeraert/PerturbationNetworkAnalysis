import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser

three2one = dict(zip(aa3, aa1))

A = nx.read_gpickle('tests/test.p')
structure = PDBParser().get_structure('X', '/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent')[0]

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = atom.serial_number

for u, v in A.edges():
    print(node2CA[u], node2CA[v])