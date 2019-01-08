from CreateNetwork import CreateNetwork
from os import path, listdir 
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
from matplotlib.markers import MarkerStyle
import itertools
import pickle as pkl
from tqdm import tqdm
base_dir = '/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/apo'
colors = itertools.cycle(('r', 'g', 'b', 'tab:orange', 'tab:purple', 'y', 'm', 'c', 'tab:pink', 'tab:olive'))

def cutoff_needeed(frame, ref_node, cutoff_min=2, cutoff_max=14, step=0.1, L_neighbors=None):
    L_cutoffs = np.arange(cutoff_min, cutoff_max, step).tolist()
    neighbors_reached = {}
    for cutoff in L_cutoffs:
        net = CreateNetwork(cutoff=cutoff).create(frame)
        for node in net.neighbors(ref_node):
            if node not in neighbors_reached:
                if L_neighbors:
                    if node in L_neighbors:
                        neighbors_reached[node] = cutoff
                else:
                    neighbors_reached[node] = cutoff
    return neighbors_reached

lines = []
n_frames = len(listdir(base_dir))
plot_dic = {}
for i, frame in enumerate(tqdm(listdir(base_dir))):
    _ = cutoff_needeed(path.join(base_dir, frame), 'C84:H')
    for _neighbor in _:
        if _neighbor not in plot_dic:
            plot_dic[_neighbor] = 15*np.ones(n_frames)
        plot_dic[_neighbor][i] = _[_neighbor]
try:
    pkl.dump(plot_dic, open('/home/agheerae/results/time_evolution_neighbors.p', 'wb'))
except FileNotFoundError:
    pass
f = plt.figure()
for _neighbor in plot_dic:
    if _neighbor in ['V51:H', 'P10:H', 'P49:H', 'G50:H', 'G52:H', 'V8:H', 'G9:H', 'G11:H']:
        color = next(colors)
        plt.plot(range(1, len(listdir(base_dir))+1), plot_dic[_neighbor], c=color)
        lines.append(mlines.Line2D([], [], color=color, label='Node '+_neighbor))

plt.xlabel('Time in frame (1 frame=1/10ns)')
plt.ylabel('Cutoff needed to see the neighbor')
plt.ylim(0,10)
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/time_evolution_neighbors_apo.pdf')
    



        








