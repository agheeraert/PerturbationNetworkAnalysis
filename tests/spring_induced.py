import networkx as nx
import numpy as np
from os.path import dirname
import matplotlib.pyplot as plt

net = nx.read_gpickle('/home/agheerae/results/avg_sims/induced/induced_5/5_5_R249:F.p')
fig = plt.figure()
colors = nx.get_edge_attributes(net, 'color').values()
weights = nx.get_edge_attributes(net, 'weight').values()
nx.draw(net, with_labels=True, font_weight='bold', edge_width=weights, edge_color=colors, node_size=100, node_color='grey', font_size=8)
plt.savefig('/home/agheerae/results/avg_sims/induced/induced_5/5_5_R249:F_spring.pdf')
