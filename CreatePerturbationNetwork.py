from CreateNetwork import AANetwork
from DrawNetwork import DrawNetwork
from multiprocessing import Pool
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import itertools
import os
from os.path import join, dirname, isdir, basename
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

num = itertools.cycle((1, 2))

class PerturbationNetwork(AANetwork):
    """Creates the Perturbation Network"""
    def __init__(self, path1, path2, out_dir, cutoff=5):
        super().__init__(cutoff=cutoff)
        # if not avg:
        #     self.net1 = self.create(path1)
        #     self.net2 = self.create(path2)
        # # elif std:
        # #     self.net1 = self.create_avg_std(path1)
        # #     self.net2 = self.create_avg_std(path2)
        # else:
        #     self.net1 = self.create_average(path1)
        #     self.net2 = self.create_average(path2)
        self.out_dir = out_dir
        pool = Pool(processes=2)
        self.net1, self.net2 = pool.map(self.smart_loader, [path1, path2])          
        # self.net1 = smart_loader(path1)
        # self.net2 = smart_loader(path2)

    def smart_loader(self, path):
        if len(basename(path)) != 0: #should get the name of a file or a folder without /
            string = basename(path).split('.')[0]+'.p'
        elif len(basename(path[:-1])) !=0: #should get the name of a folder with /
            string = basename(path[:-1])+'.p'
        else: #dummy name generation
            string = 'aa_net%s.p' %next(num)
        w_dir = join(self.out_dir, 'aa_net')
        try: 
            os.makedirs(w_dir)
        except FileExistsError:
            pass
        if isdir(path): 
            #loads a folder and do the average aanet, saves
            net = self.create_average(path)
            self.save(join(w_dir, string))
        elif path[-4:] == '.pdb': 
            #loads a single frame, do the aanet, saves
            net = self.create(path)
            self.save(join(w_dir, string))
        elif path[-2:] == '.p': 
            #directly loads a precomputed network
            net = nx.read_gpickle(path)
        else:
            print('Unknown extension for file %s' %path)
            return None
        return net
    
    def align(self, L_pos):
        assert len(L_pos)%4==0, print('Number of flags for alignment must be multiple of 4')
        apo_start_flags = L_pos[::4]
        apo_finish_flags = L_pos[2::4]
        holo_start_flags = L_pos[1::4]
        holo_finish_flags = L_pos[3::4]
        print(apo_start_flags, apo_finish_flags, holo_start_flags, holo_finish_flags)
        L_apo, L_holo = [], []
        def count(net, L, start_flags, finish_flags):
            counting = False
            for elt in net.nodes():
                if elt in start_flags:
                    counting = True
                if counting:
                    L.append(elt)
                if elt in finish_flags:
                    counting = False
        count(self.net1, L_apo, apo_start_flags, apo_finish_flags)
        count(self.net2, L_holo, holo_start_flags, holo_finish_flags)
        holo2apo = dict(zip(L_holo, L_apo))
        self.net1 = self.net1.subgraph(L_apo)
        self.net2 = nx.relabel_nodes(self.net2.subgraph(L_holo), holo2apo)
        
    def perturbation(self, threshold=0):
        """Creates the perturbation network between the initialized networks at a specific threshold"""
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
        self.perturbation = perturbation

    def draw_perturbation(self, pdb_path=None, method='default', colors=['red', 'dodgerblue']):
        """Draws and save the perturbation network for the initialized networks for all integer 
        thresholds until the network is empty"""
        self.perturbation()
        DrawNetwork(self.perturbation, self.out_dir, pdb_path, method=method, colors=['red', 'dodgerblue'])