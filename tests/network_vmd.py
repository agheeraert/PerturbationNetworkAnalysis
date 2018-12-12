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
                    help='network file')
parser.add_argument('output',  type=str,
                    help='output file')                   
parser.add_argument('-nc',  type=bool, default=False,
                    help='the network has no edge color attribute')      
parser.add_argument('--node', type=bool, default=False,
                    help="to draw node networks")              
args = parser.parse_args()

three2one = dict(zip(aa3, aa1))
three2one['5CS'] = 'C'
A = nx.read_gpickle(args.net_path)
structure = PDBParser().get_structure('X', args.pdb_path)[0]
color = not args.nc

if not args.node:
    div = max(nx.get_edge_attributes(A, 'weight').values())/1.5

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords

with open(args.output, 'w') as output:
        output.write('draw delete all \n draw color 1 \n')
        red = None
        for u, v in A.edges():
            if color:
                if red is None:
                    red = A.get_edge_data(u, v)['color']=='r'
                if A.get_edge_data(u, v)['color']=='g'and red:
                    output.write('draw color 7 \n')
                    red=False
                elif A.get_edge_data(u, v)['color']=='r'and not red:
                    output.write('draw color 1 \n')
                    red=True
            if not args.node:
                output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
            else:
                output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+'1 \n')
        for u in A.nodes():
            if args.node:
                if A.nodes(data=True)[u]['color']=='g':
                    output.write('draw color 7 \n')
                elif A.nodes(data=True)[u]['color']=='r':
                    output.write('draw color 1 \n')
            try:
                output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius 1.5 \n')
            except KeyError:
                pass