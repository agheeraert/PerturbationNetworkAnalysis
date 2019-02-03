from os.path import join as jn
import networkx as nx
import numpy as np

sims = ['/home/agheerae/results/sim'+str(i)+'/nets/' for i in range(1,5)]

node1, node2, color = input('Type the edge that interest you (node node color) with color=(r/b)')

apo_values, prfar_values = [], []

for sim in sims:
    apo_values.append(nx.read_gpickle(jn(sim, 'apo', '5.p')).get_edge_data(node1, node2))
    prfar_values.append(nx.read_gpickle(jn(sim, 'prfar', '5.p')).get_edge_data(node1, node2))

diff = np.zeros([4,4])

for i, value in apo_values:
    for j, value2 in prfar_values:
        diff[i, j] = (value2 - value1)*(color=='b')

print(np.argmax(diff))









