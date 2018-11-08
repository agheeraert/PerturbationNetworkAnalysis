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

class CreateNetwork:
    def __init__(self, pos1=None, pos2=None, cutoff=5):
        self.pos1 = pos1
        self.pos2 = pos2
        self.three2one = dict(zip(aa3, aa1))
        self.one2three = dict(zip(aa1, aa3))
        self.cutoff = cutoff

    def create(self, pdb):
        mol = bg.Pmolecule(pdb)
        self.net = mol.network(cutoff=self.cutoff, weight=True)
        self.structure = PDBParser().get_structure('X', pdb)[0]
        if self.pos1 and self.pos2:
            for node in list(self.net.nodes):
                pos = int(node[1::])
                if pos not in range(self.pos1, self.pos2):
                    self.net.remove_node(node)

        residues, bad_nodes = [], []
        for residue in self.structure.get_residues():
            if residue.resname in self.three2one:
                residues.append(self.three2one[residue.resname])
            else:
                bad_nodes.append(residue.parent.id+str(residue.id[1]))
        self.net.remove_nodes_from(bad_nodes)
        old_labels = self.net.nodes
        labels = [a+b[1:]+':'+b[0] for a,b in zip(residues, old_labels)]
        mapping = dict(zip(old_labels, labels))
        self.net = nx.relabel_nodes(self.net, mapping)
        return self.net

    def create_average(self, folder):
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
    
    def draw_avg(self, folder):
        net = self.create_average(folder)
        f = plt.figure()
        nx.draw(net, with_labels=True, font_weight='bold')
        plt.savefig('/home/hgheerae/Python/avg.pdf')    

    def draw(self, net):
        net = self.create(pdb)
        nx.draw(net, with_labels=True, font_weight='bold')
        plt.show()

if __name__ == '__main__':
    CreateNetwork().draw_avg('/home/aria/PerturbationNetworkAnalysis/data/frames_apo/')
