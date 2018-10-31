import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser

three2one = dict(zip(aa3, aa1))
one2three = dict(zip(aa1, aa3))
# path = '/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'
# path = '/home/aria/Stage4A_partie1/python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'
path = '/home/aria/PerturbationNetworkAnalysis/data/frame1000/frame_1000_apo.pdb'

A = nx.read_gpickle('tests/test.p')
structure = PDBParser().get_structure('X', path)[0]
L_residues = []
for residue in structure.get_residues():
    if residue.resname in three2one:
        L_residues.append(three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id)

for residue in L_residues:
    if residue not in A.nodes:
        A.add_node(residue)

print(len(A.nodes()))

B = nx.to_numpy_matrix(A)
B= B/np.max(B)
print(B)
np.savetxt('tests/test.dat', B, fmt='%-1f', header='name CA', comments='')