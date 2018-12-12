import biographs as bg
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
from os import listdir
from os.path import join, dirname
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

class AlignPerturbation():
    def __init__(self, pdb1, pdb2, id1, id2, alnpath, cutoff=5):
        super(AlignPerturbation, self).__init__()
        self.id1 = id1
        self.id2 = id2
        self.pdb1 = pdb1
        self.pdb2 = pdb2
        self.alnpath = alnpath
        self.cutoff = cutoff
        self.three2one = dict(zip(aa3, aa1))
        self.one2three = dict(zip(aa1, aa3))
        self.three2one['5CS'] = 'C'

    def create(self, pdb):
        mol = bg.Pmolecule(pdb)
        self.net = mol.network(cutoff=self.cutoff, weight=True)
        self.structure = PDBParser().get_structure('X', pdb)[0]
        residues, bad_nodes = [], []
        for residue in self.structure.get_residues():
            residues.append(self.three2one[residue.resname])
        old_labels = self.net.nodes
        labels = [a+b[1:]+':'+b[0] for a,b in zip(residues, old_labels)]
        mapping = dict(zip(old_labels, labels))
        self.net = nx.relabel_nodes(self.net, mapping)
        return self.net
    
    def remove_dic(self, dic, char='-'):
        toremove = []
        for elt in dic:
            if elt[0]==char or dic[elt][0]==char:
                toremove.append(elt)
        for elt in toremove:
            del dic[elt]
        return dic

    def alnParser(self):
        with open(self.alnpath, 'r') as aln_file:
            chain1, chain2, one2two = None, None, None
            for line in aln_file:
                words = line.split()
                if words:
                    if words[0:2] == ['Chain', '1:']:
                        chain1 = [aa + str(i) + ':' + self.id1 for i, aa in enumerate(list(words[3]), int(words[2]))]
                    if words[0:2] == ['Chain', '2:']:
                        chain2 = [aa + str(i) + ':' + self.id2 for i, aa in enumerate(list(words[3]), int(words[2]))]
                if chain1 and chain2:
                    if not one2two:
                        one2two = dict(zip(chain1, chain2))
                    else:
                        one2two.update(dict(zip(chain1, chain2)).items())
        one2two = self.remove_dic(one2two)
        two2one = dict(zip(one2two.values(), one2two.keys()))
        return one2two, two2one

    def align_networks(self):
        net1, net2 = self.create(self.pdb1), self.create(self.pdb2)
        one2two, two2one = self.alnParser()
        toremove = []
        for node in net1.nodes():
            if node not in two2one:
                toremove.append(node)
        net1.remove_nodes_from(toremove)
        toremove = []
        for node in net2.nodes():
            if node not in one2two:
                toremove.append(node)
        net2.remove_nodes_from(toremove)
        for elt in one2two:
            if elt not in net1.nodes():
                net1.add_node(elt) 
        for elt in two2one:
            if elt not in net2.nodes():
                net2.add_node(elt)       
        net2 = nx.relabel_nodes(net2, one2two)
        net1 = net1
        return net1, net2

    def perturbation(self, threshold):
        net1, net2 = self.align_networks()
        perturbation = nx.compose(net1, net2)
        mapping_weight, mapping_color = {}, {}
        for u, v in set(net2.edges) & set(net1.edges):
            w = net2.get_edge_data(u, v)['weight'] - net1.get_edge_data(u,v)['weight']
            if abs(w) > threshold:
                mapping_weight[(u, v)] = abs(w)
                if w > threshold:
                    mapping_color[(u, v)] = 'r'
                else:
                    mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(net2.edges) - set(net1.edges):
            weight_edge = net2.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'r'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(net1.edges) - set(net2.edges):
            weight_edge = net1.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        nx.set_edge_attributes(perturbation, name='color', values=mapping_color)
        nx.set_edge_attributes(perturbation, name='weight', values=mapping_weight)
        perturbation.remove_nodes_from(list(nx.isolates(perturbation)))
        return perturbation

    def draw_perturbation(self, threshold, output):
        if isinstance(threshold, int) or isinstance(threshold, float):
            threshold = [threshold]
        for elt in threshold:
            net = self.perturbation(elt)
            colors = nx.get_edge_attributes(net, 'color').values()
            weights = nx.get_edge_attributes(net, 'weight').values()
            fig = plt.figure()
            nx.draw(net, with_labels=True, font_weight='bold', edge_width=weights, edge_color=colors, node_size=100, node_color='grey', font_size=8)  
            if len(net.edges) > 0: 
                plt.savefig(join(output, str(elt)+'.pdf'))
                nx.write_gpickle(net, join(output, str(elt)+'.p'))


if __name__ == '__main__':
    for chain in ['C', 'D']:
        for cutoff in range(3, 10):
            AlignPerturbation(pdb1 = '/home/agheerae/Python/PerturbationNetworkAnalysis/data/aligned_tmaritima_yeast/1GPW_apo_x-ray.pdb', pdb2='/home/agheerae/Python/PerturbationNetworkAnalysis/data/aligned_tmaritima_yeast/apo_1OX6_chainA_alig.pdb', id1='A', id2=chain, alnpath='/home/agheerae/Python/PerturbationNetworkAnalysis/data/aligned_tmaritima_yeast/'+chain+'.aln', cutoff=cutoff).draw_perturbation(threshold=range(0,100), output='/home/agheerae/results/comp_bacteria_yeast/'+chain+'/cutoff_'+str(cutoff))

