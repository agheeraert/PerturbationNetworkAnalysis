from CreatePerturbationNetwork import CreatePerturbationNetwork
from CreateNetwork import CreateNetwork
from Bio.PDB.Polypeptide import aa1, aa3
import networkx as nx
import matplotlib.pyplot as plt


class AlignPerturbation(CreatePerturbationNetwork):
    def __init__(self, pdb1, pdb2, id1, id2, alnpath):
        super(AlignPerturbation, self).__init__(pdb1, pdb2)
        self.id1 = id1
        self.id2 = id2
        self.pdb1 = pdb1
        self.pdb2 = pdb2
        self.alnpath = alnpath

    def alnParser(self):
        with open(self.alnpath, 'r') as aln_file:
            chain1, chain2, one2two, two2one = None, None, None, None
            for line in aln_file:
                words = line.split()
                if words:
                    if words[0:2] == ['Chain', '1:']:
                        chain1 = [aa + str(i) + ':' + self.id1 for i, aa in enumerate(list(words[3]), int(words[2]))]
                    if words[0:2] == ['Chain', '2:']:
                        chain2 = [aa + str(i) + ':' + self.id2 for i, aa in enumerate(list(words[3]), int(words[2]))]
                if chain1 and chain2:
                    if not (one2two and two2one):
                        one2two = dict(zip(chain1, chain2))
                        two2one = dict(zip(chain2, chain1))
                    if one2two and two2one:
                        one2two.update(dict(zip(chain1, chain2)).items())
                        two2one.update(dict(zip(chain2, chain1)).items())
        return one2two, two2one

    def align_networks(self):
        net1, net2 = self.create(self.pdb1), self.create(self.pdb2)
             
        one2two, two2one = self.alnParser()
        toremove = []
        for node in net1.nodes():
            if node not in one2two:
                toremove.append(node)
        net1.remove_nodes_from(toremove)
        toremove = []
        for node in net2.nodes():
            if node not in two2one:
                toremove.append(node)
        net2.remove_nodes_from(toremove)
        for elt in one2two:
            if elt not in net1.nodes():
                net1.add_node(elt) 
        for elt in two2one:
            if elt not in net2.nodes():
                net2.add_node(elt)       
        net2 = nx.relabel_nodes(net2, two2one)
        net1 = net1
        return net1, net2

    def aligned_perturbation(self, threshold):
        self.align_networks()
        return super(AlignPerturbation, self).perturbation(threshold)
    def draw_aligned_perturbation(self, threshold, output):
        self.align_networks()
        return super(AlignPerturbation, self).draw_perturbation(threshold, output)

if __name__ == '__main__':
    print(AlignPerturbation(pdb1 = '/home/hgheerae/Python/PerturbationNetworkAnalysis/TmaritimaVSyeast/1GPW_apo_x-ray.pdb', pdb2='/home/hgheerae/Python/PerturbationNetworkAnalysis/TmaritimaVSyeast/apo_1OX6_chainA_alig.pdb', id1='H', id2='H', alnpath='/home/hgheerae/Python/PerturbationNetworkAnalysis/tests/test.aln').draw_aligned_perturbation(threshold=40, output='/home/hgheerae/Python/PerturbationNetworkAnalysis/TmaritimaVSyeast/res.pdf')) 

