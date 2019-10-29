import biographs as bg
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
from os import listdir
from os.path import join
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

three2one = dict(zip(aa3, aa1))
one2three = dict(zip(aa1, aa3))
# Some residues have different name in MD simulations (because of the protonation or other issues)
# this can lead to non recognition of some residues. Add a line here if needed.
# (See Troubleshooting 4.1 from tutorial for more information)
three2one['5CS'] = 'C'
three2one['HIP'] = 'H'
three2one['HID'] = 'H'
three2one['HIE'] = 'H'
three2one['GLH'] = 'E'
three2one['ASH'] = 'D'
three2one['S2P'] = 'S'

class CreateNetwork:
    def __init__(self, pos1=None, pos2=None, cutoff=5):
        self.pos1 = pos1
        self.pos2 = pos2
        self.three2one = three2one
        self.one2three = one2three 
        self.cutoff = cutoff

    def create(self, pdb):
        """ Creates the amino acid network using biographs"""
        mol = bg.Pmolecule(pdb)
        self.net = mol.network(cutoff=self.cutoff, weight=True)
        self.structure = PDBParser().get_structure('X', pdb)[0]
        if self.pos1 and self.pos2:
            for node in list(self.net.nodes):
                pos = int(node[1::])
                if pos not in range(self.pos1, self.pos2):
                    self.net.remove_node(node)

        residues = []
        for residue in self.structure.get_residues():
            residues.append(self.three2one[residue.resname])
        old_labels = self.net.nodes
        labels = [a+b[1:]+':'+b[0] for a,b in zip(residues, old_labels)]
        mapping = dict(zip(old_labels, labels))
        self.net = nx.relabel_nodes(self.net, mapping)
        return self.net

    def save(self, pdb, output):
        """Directtly saves the network"""
        net = self.create(pdb)
        nx.write_gpickle(net, output)

    def create_average(self, folder):
        """Creates the average network over a folder. The average is computed at each step for memory issues.
        This has been demonstrated to be equivalent to a global average."""
        net = None
        weights = dict()
        for filepath in tqdm(listdir(folder)):
            if net:
                _net = self.create(join(folder, filepath))
                for u, v in _net.edges():
                    if (u, v) in weights:
                        weights[(u,v)] += _net.get_edge_data(u, v)['weight']
                    else:
                        weights[(u,v)] = _net.get_edge_data(u, v)['weight']
                net = nx.compose(net, _net)
            else:
                net = self.create(join(folder, filepath))
                for u, v in net.edges():
                    weights[(u, v)] = net.get_edge_data(u, v)['weight']
        new_weights = {}
        for elt in weights:
            new_weights[elt]=weights[elt]/len(listdir(folder))
        nx.set_edge_attributes(net, name='weight', values=new_weights)
        return net

    def save_avg(self, folder, output):
        """Saves the average"""
        net = self.create_average(folder)
        nx.write_gpickle(net, output) 

    def draw_avg(self, folder, output):
        """Draws and saves the average"""
        net = self.create_average(folder)
        f = plt.figure()
        nx.draw(net, with_labels=True, font_weight='bold')
        nx.write_gpickle(net, output)    

    def draw(self, pdb, output):
        """Draws the network"""
        net = self.create(pdb)
        nx.draw(net, with_labels=True, font_weight='bold')
