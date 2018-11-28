import networkx as nx
from CreatePerturbationNetwork import CreatePerturbationNetwork
from os import path

ref_network_path = ''
comp_network_path = ''
output_folder = ''
threshold = range(0, 10)

for threshold in thresh_range:
    ref_net = nx.read_gpickle(ref_network_path)
    comp_net = nx.read_gpickle(comp_network_path)
    perturbation = nx.compose(ref_net, comp_net)
    mapping_color = {}

    for node in set(self.net2.nodes) & set(self.net1.nodes):
        k = net2.degree(node) - net1.degree(node)
        if k > threshold:
            mapping_color[node] = 'r'
        elif abs(k) < threshold:
            perturbation.remove_node(node)
        else:
            mapping_color[node] = 'g'

    for u, v in set(self.net2.edges) - set(self.net1.edges):
        k = net.degree(node)
        if k > threshold:
            mapping_color[node] = 'r'
        else:
            perturbation.remove_node(node)

    for u, v in set(self.net1.edges) - set(self.net2.edges):
        k = net1.degree(node)
        if k > threshold:
            mapping_color[node] = 'g'
        else:
            perturbation.remove_node(node)
    nx.set_node_attributes(perturbation, name='color', values=mapping_color)    
    nx.write_gpickle(path.join(output_folder, 'thresh_'+str(thresh)+'.p'))

    colors = nx.get_edge_attributes(net, 'color').values()
    c1, c2, positions = 0, 0, {}
    for node in net.nodes():
        if node[-1] == 'H':
            positions[node] = [c1%2, (c1//2-(c1+1)%2)]
            c1+=1
        elif node[-1] == 'F':
            positions[node] = [c2%2+3, (c2//2-(c1+1)%2)]
            c2+=1
    fig = plt.figure()
    nx.draw(net, with_labels=True, font_weight='bold', node_color=colors, pos=positions, node_size=100, font_size=8)
    plt.savefig(path.join(output_folder, 'thresh_'+str(thresh)+'.pdf'))

