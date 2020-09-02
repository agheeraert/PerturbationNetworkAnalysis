from CreateNetwork import AANetwork, AANetwork2, three2one
from DrawNetwork import DrawNetwork
from multiprocessing import Pool
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import itertools
import os
import pickle as pkl
from os.path import join, dirname, isdir, basename
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB import PDBParser
warnings.simplefilter('ignore', PDBConstructionWarning)

num = itertools.cycle((1, 2))

class PerturbationNetwork2():
    """Creates Perturbation Network between two networks"""
    def __init__(self, L_path1, L_path2, out_dir):
        self.out_dir = out_dir
        
        pool = Pool(processes=2) 
        AANet1, AANet2 = pool.map(self.smart_loader, [L_path1, L_path2])
        mat1, mat2 = AANet1.adjacency, AANet2.adjacency
        self.atom2resname = AANet1.atom2resname
        self.pertmat = np.abs(mat2 - mat1)
        self.sign = (mat2 - mat1 > 0)
        self.net = self.reconstruct_net(self.pertmat, self.sign)
    
    def reconstruct_net(self, mat, sign):
        mat = np.stack([mat, sign], axis=-1)
        dt=[('weight',mat.dtype),('color', mat.dtype)]
        mat = mat.view(dtype=dt).reshape(mat.shape[:-1])
        net = nx.from_numpy_matrix(mat)
        net = nx.relabel_nodes(net, self.atom2resname)
        net.remove_edges_from([edge for edge, attrs in dict(net.edges()).items() if attrs['weight'] == 0])
        net.remove_nodes_from(list(nx.isolates(net)))
        return net


    def smart_loader(self, L_path):
        path = L_path[0]
        if len(basename(path)) != 0: #should get the name of a file or a folder without /
            string = basename(path).split('.')[0]+'.p'
        elif len(basename(path[:-1])) !=0: #should get the name of a folder with /
            string = basename(path[:-1])+'.p'
        else: #dummy name generation
            string = 'aa_net%s.p' %next(num)
        if path[-2:] != '.p':
            try:
                net = AANetwork2(L_path[:-1], topo=L_path[-1])
            except AttributeError:
                print(L_path[:-1], 'ERRRRROR')
            net.save(join(self.out_dir, string))
        else: 
            #directly loads a precomputed network
            net = nx.read_gpickle(path)
        return net   

    def save(self):
        """Saves the PertNetwork class"""
        pkl.dump(self, open(join(self.out_dir, 'pertnet.p'), 'wb')) 

    def iterator(self, iterator=1, drawing_method='default', pdb_drawing=None, L_OXY=None, chaintop=None, norm=1.5):
        self.drawing_method = drawing_method
        self.pdb_drawing = pdb_drawing
        self.L_OXY = L_OXY
        self.chaintop = chaintop
        self.norm = norm
        
        pertmat = self.pertmat
        empty = np.sum(self.pertmat) == 0
        threshold = 0
        self.get_2dpos()
        self.get_3dpos()
        self.div = max(nx.get_edge_attributes(self.net, 'weight').values())/norm
        while not empty:
            pertmat[pertmat<threshold] = 0
            pertmat
            empty = np.sum(self.pertmat) == 0
            if not empty:
                net = self.reconstruct_net(pertmat, self.sign)
                self.draw(net, join(self.out_dir, '%s.pdf' % threshold))
                self.vmd(net, join(self.out_dir, '%s.tcl' % threshold))
            threshold+=iterator
    
    def set_drawing_parameters(self, colors=['red', 'blue'], node_size=100, font_size=10):
        self.colors = colors
        self.font_size = font_size
        self.node_size = node_size

    def get_2dpos(self):
        """ Draws the network using the O X Y method, O, X, Y are 3d coordinate points that should 
        frame the representation"""
        if self.drawing_method == 'oxy':
            O, X, Y = self.L_OXY[:3],self.L_OXY[3:6],self.L_OXY[6:]
            Ox = [c1-c2 for c1, c2 in zip(O, X)]   
            Oy = [c1-c2 for c1, c2 in zip(O, Y)]
            normOx = np.linalg.norm(Ox)  
            normOy = np.linalg.norm(Oy)
            structure = PDBParser().get_structure('X', self.pdb_drawing)[0]
            pos = {}
            distance_thresh = 1
            for atom in structure.get_atoms():
                if atom.id == 'CA':
                    residue = atom.parent
                    if self.chaintop:
                        c = 1*(residue.parent.id == self.chaintop)
                    else:
                        c = 0
                    if residue.resname in three2one:
                            "we compute the coordinates using the distance of the point to the lines of the frame"
                            AO = [c1-c2 for c1, c2 in zip(O, atom.coord)]
                            x = np.linalg.norm(np.cross(AO,Ox))/normOx
                            y = np.linalg.norm(np.cross(AO,Oy))/normOy+1*c
                            pos[three2one[residue.resname]+str(residue.id[1])] = (x, y)
            self.pos = pos
        else:
            self.pos = None
    
    def get_3dpos(self):
        structure = PDBParser().get_structure('X', self.pdb_drawing)[0]
        node2CA = {}
        for atom in structure.get_atoms():
            if atom.id == 'CA':
                residue = atom.parent
                if residue.resname in three2one:
                        coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                        node2CA[three2one[residue.resname]+str(residue.id[1])] = coords
        self.node2CA = node2CA

    def draw(self, net, output):
        fig = plt.figure()
        colors = list(nx.get_edge_attributes(net, 'color').values())
        for i, color in enumerate(colors):
            if color == 1:
                colors[i] = self.colors[0]
            if color == 0:
                colors[i] = self.colors[1]
        width = list(nx.get_edge_attributes(net, 'weight').values())
        max_width = max(width)
        for i, elt in enumerate(width):
            width[i] = elt/max_width*5
        # weights = {(u, v): round(nx.get_edge_attributes(net, 'weight')[(u,v)]) for (u, v) in nx.get_edge_attributes(net, 'weight')}
        if self.pos:
            pos = {node: self.pos[node] for node in net.nodes()}
            nx.draw(net, font_weight='bold', labels={node: node.split(':')[0] for node in net.nodes()}, width=width, pos=pos, edge_color=colors, node_size=self.node_size, node_shape='o', font_size=self.font_size, node_color='lightgrey')
        else:
            nx.draw(net, font_weight='bold', labels={node: node.split(':')[0] for node in net.nodes()}, width=width, edge_color=colors, node_size=self.node_size, node_shape='o', font_size=self.font_size, node_color='lightgrey')
        plt.savefig(output)
        plt.close()
    
    def vmd(self, net, output):
        with open(output, 'w') as output:
            output.write('draw delete all \n')
            previous = None
            for u, v in net.edges():
                c = net.get_edge_data(u, v)['color']
                if c:
                    c = 'blue'
                else:
                    c = 'red' 
                if previous != c:
                    output.write('draw color ' +c+'\n')
                    previous = c  
                try:
                    output.write('draw cylinder { ' + str(self.node2CA[u]) + ' '+ ' } ' + '{ ' + str(self.node2CA[v]) + ' '+ ' } radius '+str(net.get_edge_data(u, v)['weight']/self.div)+' \n')
                except KeyError:
                    pass
            output.write('draw color silver \n')
            for u in net.nodes():
                try:
                    output.write('draw sphere { ' + str(self.node2CA[u]) + ' '+ ' } radius '+str(self.norm)+' \n')
                except KeyError:
                    print('Warning, residue', u, 'probably mutated between the two networks')


    
            
            

if __name__ == '__main__':
    # PerturbationNetwork2(['/home/aria/landslide/MDRUNS/YEAST/all_trajs/model1_apo.dcd', '/home/aria/tmp/topo.pdb'], ['/home/aria/landslide/MDRUNS/YEAST/all_trajs/model1_holo_protein.dcd', '/home/aria/tmp/topo.pdb'], 'results/test_newpn/')   
    # pn2 = PerturbationNetwork2(['results/test_newpn/model1_apo.p'], ['results/test_newpn/model1_holo_protein.p'], 'results/test_newpn/')   
    pn2 = PerturbationNetwork2(['/home/aria/landslide/MDRUNS/YEAST/all_trajs/model%s_apo.dcd' %i for i in ['1', 'a']+list(range(3,7))]+['/home/aria/tmp/topo.pdb'], ['/home/aria/landslide/MDRUNS/YEAST/all_trajs/model%s_holo_protein.dcd' %i for i in ['1', 'a']+list(range(3,7))]+['/home/aria/tmp/topo.pdb'], '/home/aria/landslide/RESULTS/YEAST/all/')   

    pn2.save()
    pn2.set_drawing_parameters()
    pn2.iterator(drawing_method='oxy', pdb_drawing='/home/aria/tmp/topo.pdb', L_OXY=[44, 16, 57.5, 50.6, 75.5, 75.7, 32.26, 61.28, 14.43])

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
    
    def align_prot_mol(self, chainname, molname):
        """Aligns coarsely a molecule to a protein fragment
        carful the first network should be the one with the protein
        and the second one with the molecule. use the align function to align molecules"""
        firstu = False
        for v in self.net1.nodes():
            if v[-1]==chainname:
                if not firstu:
                    firstu = True
                    u = v
                else:
                    self.net1 = nx.contracted_nodes(self.net1, u, v)
        self.net1 = nx.relabel_nodes(self.net1, {firstu: chainname})


        
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