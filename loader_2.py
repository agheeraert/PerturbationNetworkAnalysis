import networkx as nx
from CreatePerturbationNetwork import CreatePerturbationNetwork
from os import path
import argparse
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Constructs node degree perturbation network from 0-threshold edge weight perturbation network')
parser.add_argument('path1',  type=str,
                    help='Static input file path of reference (apo)')
parser.add_argument('path2',  type=str,
                    help='Static input file path to compare (complex)')
parser.add_argument('output',  type=str,
                    help='Output directory for the results')
parser.add_argument('range',  type=float, nargs=3,
                    help='Range within which draw the perturbation networks')
args = parser.parse_args()


ref_network_path = args.path1
comp_network_path = args.path2
output_folder = args.output
thresh_range = np.arange(args.range[0], args.range[1], args.range[2]).tolist()

for threshold in thresh_range:
    try:
        str_thresh = str(int(threshold))
    except ValueError:
        str_thresh = str(threshold).replace('.', '-')

    ref_net = nx.read_gpickle(ref_network_path)
    comp_net = nx.read_gpickle(comp_network_path)
    perturbation = nx.compose(ref_net, comp_net)
    mapping_color = {}

    for node in comp_net.nodes:
        k = comp_net.degree(node) - ref_net.degree(node)
        if abs(k) >= threshold:
            if k >= 0:
                mapping_color[node] = 'r'
            else:
                mapping_color[node] = 'g'
        else:
            perturbation.remove_node(node)

    nx.set_node_attributes(perturbation, name='color', values=mapping_color) 

    if len(perturbation.nodes) != 0:   
        nx.write_gpickle(perturbation, path.join(output_folder, 'thresh_'+str_thresh+'.p'))
        colors = list(nx.get_node_attributes(perturbation, 'color').values())
        c1, c2, positions = 0, 0, {}
        for node in perturbation.nodes():
            if node[-1] == 'H':
                positions[node] = [c1%2, (c1//2-(c1+1)%2)]
                c1+=1
            elif node[-1] == 'F':
                positions[node] = [c2%2+3, (c2//2-(c1+1)%2)]
                c2+=1
        fig = plt.figure()
        nx.draw(perturbation, with_labels=True, font_weight='bold', pos=positions, node_color=colors, node_size=100, font_size=8)
        plt.savefig(path.join(output_folder, 'thresh_'+str_thresh+'.pdf'))

