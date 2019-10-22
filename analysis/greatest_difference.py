from os.path import join as jn
import networkx as nx
import numpy as np

sims = ['/home/agheerae/results/sim'+str(i)+'/nets/' for i in range(1,5)]

node1, node2, color = input('Type the edge that interest you (node;node;color) with color=(r/b)').split(';')

apo_values, prfar_values = [], []

for sim in sims:
    try:
        apo_values.append(nx.read_gpickle(jn(sim, 'apo', '5.p')).get_edge_data(node1, node2)['weight'])
    except:
        apo_values.append(0)
    try:
        prfar_values.append(nx.read_gpickle(jn(sim, 'prfar', '5.p')).get_edge_data(node1, node2)['weight'])
    except:
        prfar_values.append(0)
diff = np.zeros([4,4])
print(apo_values, prfar_values)
for i, value1 in enumerate(apo_values):
    for j, value2 in enumerate(prfar_values):
        diff[i, j] = (value1 - value2)*(color=='b')

print(diff)
