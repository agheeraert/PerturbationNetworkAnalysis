import networkx as nx
import argparse


parser = argparse.ArgumentParser(description='Remap the nodes of the graph with H=>F and vice versa')
parser.add_argument('path',  type=str,
                     help='file to permute')
args = parser.parse_args()

net = nx.read_gpickle(args.path)
if 'K19:H' in net.nodes and 'C84:F' in net.nodes:
    print('Bad file...')
    mapping = {}
    for node in net.nodes():
        if node[-1] == 'H':
            mapping[node] = node[:-1]+'F'
        elif node[-1] == 'F':
            mapping[node] = node[:-1]+'H'
    net = nx.relabel_nodes(net, mapping)

if 'K19:F' in net.nodes and 'C84:H' in net.nodes:
    print('File corrected')
    nx.write_gpickle(net, args.path)