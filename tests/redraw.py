import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
import itertools
import pickle as pkl

colors = itertools.cycle(('r', 'g', 'b', 'tab:orange', 'tab:purple', 'y', 'm', 'c', 'tab:pink', 'tab:olive'))
plot_dic = pkl.load(open('/home/agheerae/results/time_evolution_neighbors.p', 'rb'))
lines = []
f = plt.figure()
for _neighbor in plot_dic:
    if _neighbor in ['V51:H', 'P10:H', 'P49:H', 'G50:H', 'G52:H', 'V8:H', 'G9:H', 'G11:H']:
        color = next(colors)
        plt.plot(range(1, 1001), plot_dic[_neighbor], c=color)
        lines.append(mlines.Line2D([], [], color=color, label='Node '+_neighbor))

plt.xlabel('Time in frame (1 frame=1/10ns)')
plt.ylabel('Cutoff needed to see the neighbor')
plt.ylim(0,14)
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/time_evolution_neighbors_2.pdf')