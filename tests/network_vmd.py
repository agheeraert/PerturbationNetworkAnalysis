import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import argparse

parser = argparse.ArgumentParser(description='file paths')
parser.add_argument('pdb_path',  type=str,
                    help='pdb file')
parser.add_argument('net_path',  type=str,
                    help='output file')
parser.add_argument('output',  type=str,
                    help='output file')
args = parser.parse_args()

three2one = dict(zip(aa3, aa1))
# path = '/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'
# path = '/home/aria/Stage4A_partie1/python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent'

A = nx.read_gpickle(args.net_path)
structure = PDBParser().get_structure('X', args.pdb_path)[0]

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords
with open(args.output, 'w') as output:
        output.write('draw delete all \n draw color 1 \n')
        for u, v in A.edges():
            if u[-1] != v[-1]:
                output.write('draw color 7 \n')
            try:
                output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/8)+' \n')
            except KeyError:
                pass
            if u[-1] != v[-1]:
                output.write('draw color 1 \n')
                
        for u in A.nodes():
            try:
                output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius 1.5 \n')
            except KeyError:
                pass