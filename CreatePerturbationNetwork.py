from CreateNetwork import CreateNetwork
import networkx as nx
import matplotlib.pyplot as plt
import os
from os.path import join, dirname
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

class CreatePerturbationNetwork(CreateNetwork):
    def __init__(self, path1, path2, pos1=None, pos2=None, avg=False):
        super(CreatePerturbationNetwork, self).__init__()
        if not avg:
            self.net1 = self.create(path1)
            self.net2 = self.create(path2)
        else:
            self.net1 = self.create_average(path1)
            self.net2 = self.create_average(path2)            

    def perturbation(self, threshold):
        perturbation = nx.compose(self.net1, self.net2)
        mapping_weight, mapping_color = {}, {}
        for u, v in set(self.net2.edges) & set(self.net1.edges):
            w = self.net2.get_edge_data(u, v)['weight'] - self.net1.get_edge_data(u,v)['weight']
            if abs(w) > threshold:
                mapping_weight[(u, v)] = abs(w)
                if w > threshold:
                    mapping_color[(u, v)] = 'r'
                else:
                    mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(self.net2.edges) - set(self.net1.edges):
            weight_edge = self.net2.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'r'
            else:
                perturbation.remove_edge(u,v)
        for u, v in set(self.net1.edges) - set(self.net2.edges):
            weight_edge = self.net1.get_edge_data(u, v)['weight']
            if weight_edge > threshold:
                mapping_weight[(u, v)] = weight_edge
                mapping_color[(u, v)] = 'g'
            else:
                perturbation.remove_edge(u,v)
        nx.set_edge_attributes(perturbation, name='color', values=mapping_color)
        nx.set_edge_attributes(perturbation, name='weight', values=mapping_weight)
        perturbation.remove_nodes_from(list(nx.isolates(perturbation)))
        return perturbation

    def draw_perturbation(self, threshold, output, rearrange=None, neighbors=True, chains=True, save=None):
        if isinstance(threshold, int) or isinstance(threshold, float):
            threshold = [threshold]
        for elt in threshold:
            net = self.perturbation(elt)
            colors = nx.get_edge_attributes(net, 'color').values()
            weights = nx.get_edge_attributes(net, 'weight').values()
            ch, cf, positions = 0, 0, {}
            if rearrange:
                for node in net.nodes():
                    if node[-1] == rearrange[0]:
                        positions[node] = [ch%2, 20*(ch//2+(ch+1)%2)]
                        ch+=1
                    elif node[-1] == rearrange[1]:
                        positions[node] = [cf%2+10, 20*(cf//2+(ch+1)%2)]
                        cf+=1
            fig = plt.figure()
            try:
                nx.draw(net, with_labels=True, font_weight='bold', edge_color=colors, pos=positions)
            except nx.NetworkXError:
                nx.draw(net, with_labels=True, font_weight='bold', edge_color=colors)   
            plt.savefig(output+'_'+str(elt)+'.pdf')
            output_folder = dirname(output)
            if isinstance(save, str):
                nx.write_gpickle(net, join(output_folder, save+'_'+str(elt)+'.p'))        
            if neighbors:
                if 'neighbors_'+str(elt) not in os.listdir(output_folder):    
                    os.mkdir(join(output_folder, 'neighbors_'+str(elt)))
                if chains:
                    if '4D_'+str(elt) not in os.listdir(output_folder):    
                        os.mkdir(join(output_folder, '4D_'+str(elt)))
                for node in net.nodes():
                    subG = nx.Graph()
                    subG.add_nodes_from(net.neighbors(node))
                    fourD = False
                    for _node in net.neighbors(node):
                        subG.add_edge(node, _node)
                        if _node[-1] != node[-1]:
                            fourD = True
                    nx.set_edge_attributes(subG, name='color', values=nx.get_edge_attributes(net, 'color')) 
                    nx.set_edge_attributes(subG, name='weight', values=nx.get_edge_attributes(net, 'weight')) 
                    _colors = nx.get_edge_attributes(subG, 'color').values()
                    _weights = nx.get_edge_attributes(subG, 'weight').values()
                    fig = plt.figure()
                    nx.draw(subG, with_labels=True, font_weight='bold', edge_color=_colors, width=[w / 5 for w in _weights])
                    plt.savefig(join(output_folder, 'neighbors_'+str(elt), node+'.pdf'))
                    if fourD and chains:
                        plt.savefig(join(output_folder, '4D_'+str(elt), node+'.pdf'))


if __name__ == '__main__':
    # CreatePerturbationNetwork(path1='/home/hgheerae/Python/script_lorenza_selection/ttr_v4_3djz_L55P/pdb/pdb3djz.ent', path2='/home/hgheerae/Python/script_lorenza_selection/ttr_v4_1f41_WT_HUM/pdb/pdb1f41.ent',pos1=11, pos2=125).draw_perturbation(threshold=4, output='/home/hgheerae/Python/script_lorenza_selection/res.pdf', rearrange=('A', 'B'), save='test.p')
    # CreatePerturbationNetwork(path1='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/frames_apo/', path2='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/frames_prfar/', avg=True).draw_perturbation(threshold=40, output='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/avg.pdf', rearrange=('H', 'F'), save='avg')
    CreatePerturbationNetwork(path1='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/frames_apo_noH/', path2='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/frames_prfar_noH', avg=True).draw_perturbation(threshold=range(1, 30, 1), output='/home/hgheerae/Python/PerturbationNetworkAnalysis/data/avg/', rearrange=('H', 'F'), save='avg')
