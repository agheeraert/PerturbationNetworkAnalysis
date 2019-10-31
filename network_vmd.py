import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import argparse

from CreateNetwork import three2one

parser = argparse.ArgumentParser(description='file paths')
parser.add_argument('pdb_path',  type=str,
                    help='pdb file')
parser.add_argument('f',  type=str, nargs='+',
                    help='network file. Can be a list but they should all relate to the same pdb structure.')
parser.add_argument('-nc',  type=bool, default=False,
                    help='to specify that the network has no edge color attribute')      
parser.add_argument('-ntodraw', type=str, nargs='+', default=None,
                    help="draw only the specified list of nodes")              
parser.add_argument('-norm', type=float, default=1.5,
                    help="changes the normalizaton factor for the edges width.")
parser.add_argument('-c', type=str, nargs=3, default=['red', 'blue', 'silver'],
                    help="changes the colors used to draw the network. In the order: edges that increases in contact number, edges that decreases in contact number, node")              

args = parser.parse_args()


structure = PDBParser().get_structure('X', args.pdb_path)[0]
color = not args.nc

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords

for fichier in args.f:
    out_path = fichier[:-2]+'.tcl'
    with open(out_path, 'w') as output:
        output.write('draw delete all \n')
        previous = None
        A = nx.read_gpickle(fichier)
        div = max(nx.get_edge_attributes(A, 'weight').values())/args.norm
        for u, v in A.edges():
            if color:
                c = A.get_edge_data(u, v)['color']
                if c == 'g':
                    c = args.c[1]
                if c == 'r':
                    c = args.c[0] 
                if previous != c:
                    output.write('draw color ' +c+'\n')
                    previous = c  
                if args.ntodraw:
                    if u in args.ntodraw and v in args.ntodraw:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                else:
                    output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
            else:
                output.write('draw color '+args.c[2]+'\n')
                if args.ntodraw:
                    if u in args.ntodraw and v in args.ntodraw:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                else:
                    output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
        output.write('draw color silver \n')
        for u in A.nodes():
            if args.ntodraw:
                if u in args.ntodraw:
                    output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(args.norm)+' \n')
            else:
                output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(args.norm)+' \n')


