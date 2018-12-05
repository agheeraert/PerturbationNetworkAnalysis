from os import path
import networkx as nx 
import itertools
import matplotlib.lines as mlines
import matplotlib.pyplot as plt 
static_folder = '/home/agheerae/results/static/'
dynamic_folder = '/home/agheerae/results/sim1/pertnet/'
markers = itertools.cycle((',', '+', '.', 'x', '*', 's', 'v')) 
L_cutoffs = range(3, 10)
L_thresh = list(range(0, 30))

def avg(filepath):
    try:
        net = nx.read_gpickle(filepath)
        weights = nx.get_edge_attributes(net, 'weight').values()
        somme = 0
        for weight in weights:
            somme+=weight
        return somme/len(weights)
    except FileNotFoundError:
        return 0


f = plt.figure()
lines = []
for cutoff in L_cutoffs:
    marker = next(markers)
    static_y, dynamic_y = [], []
    for thresh in L_thresh:
        static_y.append(avg(path.join(static_folder, 'cutoff_'+str(cutoff), str(cutoff)+'_'+str(thresh)+'.p')))
        dynamic_y.append(avg(path.join(dynamic_folder, 'cutoff_'+str(cutoff), str(cutoff)+'_'+str(thresh)+'.p')))
    plt.plot(L_thresh, static_y, c = 'g', marker=marker)
    plt.plot(L_thresh, dynamic_y, c = 'b', marker=marker)
    lines.append(mlines.Line2D([], [], color='black', marker=marker, label='Cutoff '+str(cutoff)+'Å'))

lines.append(mlines.Line2D([], [], color='g', label='Static'))
lines.append(mlines.Line2D([], [], color='b', label='Dynamic'))

plt.xlabel('Threshold')
plt.ylabel('Average weight link')
plt.legend(handles=lines)
plt.savefig('/home/agheerae/results/comp_avg_weight_per_link.pdf')