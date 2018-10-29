import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser

three2one = dict(zip(aa3, aa1))
one2three = dict(zip(aa1, aa3))

A = nx.read_gpickle('tests/test.p')
structure = PDBParser().get_structure('X', '/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent')[0]
L_residues = []
for residue in structure.get_residues():
    if residue.resname in three2one:
        L_residues.append(three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id)

for residue in L_residues:
    if residue not in A.nodes:
        A.add_node(residue)

for i in range(244-len(A.nodes())):
    A.add_node('X:'+str(i))

print(len(A.nodes()))

B = nx.to_numpy_matrix(A, dtype=int)
np.savetxt('tests/test.dat', B, fmt='%0i', header='name CA', comments='')