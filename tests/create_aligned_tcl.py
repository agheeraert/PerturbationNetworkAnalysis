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
parser.add_argument('net_path1',  type=str,
                    help='network1 file')
parser.add_argument('net_path2',  type=str,
                    help='network2 file')
parser.add_argument('output',  type=str,
                    help='output file')                   

args = parser.parse_args()

three2one = dict(zip(aa3, aa1))
three2one['5CS'] = 'C'
net1 = nx.read_gpickle(args.net_path1)
net2 = nx.read_gpickle(args.net_path2)
structure = PDBParser().get_structure('X', args.pdb_path)[0]

div = max(max(zip(nx.get_edge_attributes(net1, 'weight').values(), nx.get_edge_attributes(net2, 'weight').values())))/1.5

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords

with open(args.output, 'w') as output:
        output.write('draw delete all \n')
        for A in [net1, net2]:
            for u, v in A.edges():
                if A.get_edge_data(u, v)['color']=='g':
                        output.write('draw color 7 \n')
                elif A.get_edge_data(u, v)['color']=='r':
                        output.write('draw color 1 \n')
                try:
                    output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                except:
                    print(u, v)
        output.write('draw color 2 \n')
        for net in [net1, net2]:
            for u in net.nodes():
                output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius 1.5 \n')
