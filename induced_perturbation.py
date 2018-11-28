import networkx as nx 
from os import path, listdir
import matplotlib.pyplot as plt

perturbation_folder = '/home/agheerae/results/all_simu/cutoff/'
output_folder = '/home/agheerae/results/all_simu/induced/'
log = path.join(output_folder, 'not_in_net.log')
root_nodes = ['K19:F', 'C84:H']
L_cutoffs = list(range(3, 10))
chain1='F'
chain2='H'

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
                    nx.write_gpickle(tree, path.join(output_folder, 'induced_'+str(cutoff), _file[:-2]+'_'+root+'.p'))
                    weights = list(nx.get_edge_attributes(tree, 'weight').values())
                    colors = list(nx.get_edge_attributes(tree, 'color').values())
                    if len(weights) == 1:
                        weights = weights[0]
                        colors = colors[0]
                    c1, c2, positions = 0, 0, {}
                    for node in net.nodes():
                        if node[-1] == chain1:
                            positions[node] = [c1%2, (c1//2-(c1+1)%2)]
                            c1+=1
                        elif node[-1] == chain2:
                            positions[node] = [c2%2+3, (c2//2-(c1+1)%2)]
                            c2+=1
                    f = plt.figure()                    
                    nx.draw(tree, with_labels=True, font_weight='bold', pos=positions, edge_width=weights, edge_color=colors, node_size=100, node_color='grey', font_size=8)
                    plt.savefig(path.join(output_folder, 'induced_'+str(cutoff), _file[:-2]+'_'+root+'.pdf'))
                else:
                    with open(log, 'a') as logfile:
                        logfile.write('cutoff: '+str(cutoff)+' file: '+_file+' node: '+ root+'\n')