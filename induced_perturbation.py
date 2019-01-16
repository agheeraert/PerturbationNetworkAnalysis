import networkx as nx 
from os import path, listdir
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

three2one = dict(zip(aa3, aa1))
perturbation_folder = '/home/agheerae/results/H/pertnet/'
output_folder = '/home/agheerae/results/H/induced/'
root_nodes = ['R95:F', 'R18:H', 'K19:F', 'S225:F']
L_cutoffs = list(range(5, 6))
chain1='F'
chain2='H'
pdb = '/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/1frame_1.pdb'
structure = PDBParser().get_structure('X', pdb)[0]
pos = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                y = atom.coord[1]
                x = - sqrt(2)*atom.coord[2] + sqrt(2)*atom.coord[0]
                pos[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = (x, y)

for cutoff in L_cutoffs:
    cutoff_dir = path.join(perturbation_folder, 'cutoff_'+str(cutoff))
    for _file in listdir(cutoff_dir):
        if _file[-2:] == '.p':
            net = nx.read_gpickle(path.join(cutoff_dir, _file))
            _weights = nx.get_edge_attributes(net, 'weight')
            _colors = nx.get_edge_attributes(net, 'color')
            for root in root_nodes:
                if root in net.nodes():
                    tree = nx.Graph(nx.bfs_tree(net, root))
                    for u, v in net.edges():
                        if u in tree.nodes() and v in tree.nodes() and not (u, v) in tree.edges() and not (v, u) in tree.edges():
                            tree.add_edge(u, v)
                    nx.set_edge_attributes(tree, name='weight', values=_weights)
                    nx.set_edge_attributes(tree, name='color', values=_colors)
                    out_path = path.join(output_folder, 'induced_'+str(cutoff), _file[:-2]+'_'+root)
                    if len(tree.nodes()) !=0:
                        nx.write_gpickle(tree, out_path+'.p')
                    width = list(nx.get_edge_attributes(tree, 'weight').values())
                    colors = list(nx.get_edge_attributes(tree, 'color').values())
                    if len(colors) == 1:
                        width = [width[0]]
                        colors = [colors[0]]
                    fig = plt.figure()
                    max_width = max(width)
                    for i, elt in enumerate(width):
                        width[i] = elt/max_width*2
                    weights = {(u, v): round(nx.get_edge_attributes(tree, 'weight')[(u,v)], 2) for (u, v) in nx.get_edge_attributes(tree, 'weight')}
                    nx.draw(tree, with_labels=True, font_weight='bold', width=width, pos=pos, edge_color=colors, node_size=100, node_color='grey', font_size=8)
                    nx.draw_networkx_edge_labels(tree, pos=pos, edge_labels=weights, font_color='black', font_size=5)
                    plt.savefig(out_path+'.svg')