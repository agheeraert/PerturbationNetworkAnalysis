import networkx as nx
import numpy as np
from os.path import dirname
import matplotlib.pyplot as plt

load_path = '/home/agheerae/results/cutoff/cutoff_9/9_0.p'
new_thresh = np.arange(1, 35, 1).tolist()
output_str1 = '/home/agheerae/results/cutoff/cutoff_9/9_'
chain1 = 'F'
chain2 = 'H'


net = nx.read_gpickle(load_path)
output_folder = dirname(load_path)

for elt in new_thresh:
    to_remove_edges = []
    for u, v in net.edges():
        weight = net.get_edge_data(u, v)['weight']
        if weight < elt:
            to_remove_edges.append((u, v))
    net.remove_edges_from(to_remove_edges)
    net.remove_nodes_from(list(nx.isolates(net)))
    c1, c2, positions = 0, 0, {}
    for node in net.nodes():
        if node[-1] == chain1:
            positions[node] = [c1%2, (c1//2-(c1+1)%2)]
            c1+=1
        elif node[-1] == chain2:
            positions[node] = [c2%2+3, (c2//2-(c1+1)%2)]
            c2+=1
    fig = plt.figure()
    colors = nx.get_edge_attributes(net, 'color').values()
    weights = nx.get_edge_attributes(net, 'weight').values()
    nx.draw(net, with_labels=True, font_weight='bold', edge_width=weights, edge_color=colors, pos=positions, node_size=100, node_color='grey', font_size=8)
    thresh = str(round(elt, 1)).replace('.', '-')
    plt.savefig(output_str1+thresh+'.pdf')
    nx.write_gpickle(net, output_str1+thresh+'.p')



