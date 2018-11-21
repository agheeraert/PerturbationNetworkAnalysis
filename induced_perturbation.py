import networkx as nx 
from os import path, listdir

perturbation_folder = '/home/agheerae/results/cutoff/'
output_folder = '/home/agheerae/results/induced/'
root_nodes = ['C84:H', 'K19:F', 'H178:H', 'E180:H']
L_cutoffs = list(range(3, 10))

for cutoff in L_cutoffs:
    cutoff_dir = path.join(perturbation_folder, 'cutoff_'+str(cutoff))
    for _file in listdir(cutoff_dir):
        if _file[-2:] == '.p':
            net = nx.read_gpickle(path.join(cutoff_dir, _file))
            weights = nx.get_edge_attributes(net, 'weight')
            colors = nx.get_edge_attributes(net, 'color')
            for root in root_nodes:
                tree = nx.bfs_tree(net, root)
                for u, v in net.edges():
                    if u in tree.nodes() and v in tree.nodes() and not (u, v) in tree.edges():
                        tree.add_edge(u, v)
                nx.set_edge_attributes(tree, 'weight', weights)
                nx.set_edge_attributes(tree, 'colors', colors)
                nx.write_gpickle(path.join(output_folder, 'induced_'+str(cutoff), _file[:-2]+root))