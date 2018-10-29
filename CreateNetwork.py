import biographs as bg
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import networkx as nx
import matplotlib.pyplot as plt
from os import listdir
from os.path import join

class CreateNetwork:
    def __init__(self, pos1=None, pos2=None):
        self.pos1 = pos1
        self.pos2 = pos2
        self.three2one = dict(zip(aa3, aa1))
        self.one2three = dict(zip(aa1, aa3))

    def create(self, pdb):
        mol = bg.Pmolecule(pdb)
        self.net = mol.network(weight=True)
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
        for i, filepath in enumerate(listdir(folder)):
            if net:
                _net = self.create(join(folder, filepath))
                output_net = nx.compose(net, _net)
                new_weights = dict()
                for u, v in set(net.edges) & set(_net.edges):
                    new_weights[(u,v)] = (i-1)/float(i)*float(net.get_edge_data(u, v)['weight'])+1./i*float(_net.get_edge_data(u, v)['weight'])
                for u, v in set(net.edges) - set(_net.edges):
                    new_weights[(u,v)] = (i-1)/float(i)*float(net.get_edge_data(u, v)['weight'])
                for u, v in set(_net.edges) - set(_net.edges):
                    new_weights[(u,v)] = 1./i*float(_net.get_edge_data(u, v)['weight'])
                nx.set_edge_attributes(output_net, name='weight', values=new_weights)
                net = output_net
            else:
                net = self.create(join(folder, filepath))
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
    CreateNetwork().draw_avg('/home/hgheerae/Python/PerturbationNetworkAnalysis/data/frames_apo/')