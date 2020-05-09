import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

from CreateNetwork import three2one

class VMDNetwork:
    def __init__(self, pdb_path, nc):
        structure = PDBParser().get_structure('X', pdb_path)[0]
        self.color = not nc
        self.node2CA = {}
        for atom in structure.get_atoms():
            if atom.id == 'CA':
                residue = atom.parent
                if residue.resname in three2one:
                        coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                        self.node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords

    def create_tcl(self, net, out_path, nc=False, ntodraw=None, norm=1.5):
        with open(out_path, 'w') as output:
            output.write('draw delete all \n')
            previous = None
            div = max(nx.get_edge_attributes(net, 'weight').values())/norm
            for u, v in net.edges():
                if color:
                    c = self.net.get_edge_data(u, v)['color']
                    if c == 'g':
                        c = c[1]
                    if c == 'r':
                        c = c[0] 
                    if previous != c:
                        output.write('draw color ' +c+'\n')
                        previous = c  
                    if ntodraw:
                        if u in ntodraw and v in ntodraw:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                    else:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                else:
                    output.write('draw color '+c[2]+'\n')
                    if ntodraw:
                        if u in ntodraw and v in ntodraw:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                    else:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
            output.write('draw color silver \n')
            for u in self.net.nodes():
                if ntodraw:
                    if u in ntodraw:
                        output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(norm)+' \n')
                else:
                    output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(norm)+' \n')