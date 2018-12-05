from CreateNetwork import CreateNetwork
from os import path, listdir 
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
from matplotlib.markers import MarkerStyle
import itertools
from tqdm import tqdm
base_dir = '/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/prfar'
markers = itertools.cycle(list(MarkerStyle.markers.keys())) 

def cutoff_needeed(frame, ref_node, cutoff_min=3, cutoff_max=9, step=1, L_neighbors=None):
    L_cutoffs = range(cutoff_min, cutoff_max, step)
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

f = plt.figure()
for _neighbor in plot_dic:
    marker = next(markers)
    plt.plot(list(range(1, len(listdir(base_dir)+1))), plot_dic[_neighbor], c='b', marker=marker)
    lines.append(mlines.Line2D([], [], color='b', marker=marker, label='Node '+_neighbor))

plt.xlabel('Time in frame (1 frame=1/10ns)')
plt.ylabel('Cutoff needed to see the neighbor')
plt.ylim(0,10)
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/time_evolution_neighbors.pdf')
    



        








