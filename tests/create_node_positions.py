from Bio.PDB import PDBParser
from math import sqrt, copysign
from tqdm import tqdm
import warnings
import networkx as nx
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))

structure = PDBParser().get_structure('X', '/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/1frame_1.pdb')[0]
net = nx.read_gpickle('/home/agheerae/results/H/pertnet/cutoff_5/5_20.p')
pos = {}
L_pos = []
distance_thresh = 1.5
for atom in tqdm(structure.get_atoms()):
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                y = atom.coord[1]
                x = sqrt(2)*atom.coord[2] - sqrt(2)*atom.coord[0]
                not_placed = True
                while not_placed:
                    not_placed = False
                    for position in L_pos:
                        if abs(x-position[0]) < distance_thresh and abs(y-position[1]) < distance_thresh:
                            not_placed = True
                            x += copysign(distance_thresh, x-position[0])
                            y += copysign(distance_thresh, y-position[0])

                pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)
                L_pos.append((x, y))

nx.draw(net, with_labels=True, font_weight='bold', pos=pos, node_size=100, node_color='grey', font_size=8)
plt.show()