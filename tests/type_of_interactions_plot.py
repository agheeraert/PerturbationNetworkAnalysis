from os import listdir
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines

BASE_FOLDER = '/home/agheerae/results/avg_sims/pertnet/'
colors = itertools.cycle(('r', 'g', 'b','olive', 'purple', 'navy'))  
tuple_residues = [(['R', 'H', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y'], 'p'),
                (['G', 'A', 'P'], '?'),
                (['V', 'I', 'L', 'M', 'F', 'W', 'C'], 'h')]
dict_residues = {}
for tup in tuple_residues:
    for elt in tup[0]:
        dict_residues[elt] = tup[1]

dict_interactions = {('+', '-'): 'sb',
                    ('+', '+'): 'un_sb',
                    ('-', '-'): 'un_sb',
                    ('h', 'h'): 'h',
                    ('p', 'p'): 'p'
}

for cutoff in range(3, 10):
    folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
    if cutoff == 3:
        L_cutoffs = range(0, 11)
    else:
        L_cutoffs = range(0, 30)
    interactions_count = {'sb': ((L_cutoffs[-1]+1)*[0], "Salt bridge"),
        'un_sb': ((L_cutoffs[-1]+1)*[0], "Unknown salt bridge"),
        'h': ((L_cutoffs[-1]+1)*[0], "Hydrophobic interaction"),
        'p': ((L_cutoffs[-1]+1)*[0], "Polar-polar interaction"),
        '?': ((L_cutoffs[-1]+1)*[0], "Unknown interaction"),
        'c': ((L_cutoffs[-1]+1)*[0], "Residues covalently bound"),
                }
    c_tot = (L_cutoffs[-1]+1)*[0]
    stop = False
    for threshold in L_cutoffs:
        if cutoff == 3:
            if threshold in [0, 10]:
                string = str(threshold)[0]
            else:
                string = '0-'+str(threshold)
        else:
            string = str(threshold)
        if str(cutoff)+'_'+string+'.p' in listdir(folder):
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+string+'.p'))
            for u, v in net.edges():
                if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) == 1:
                    interaction = 'c'
                else:
                    residue_residue = (dict_residues[u[0]], dict_residues[v[0]])
                    if residue_residue in dict_interactions:
                        interaction = dict_interactions[residue_residue]
                    elif residue_residue[::-1] in dict_interactions:
                        interaction = dict_interactions[residue_residue[::-1]]
                    else:
                        interaction = '?'
                interactions_count[interaction][0][threshold]+=1
                c_tot[threshold]+=1
        elif not stop:
            stopping_point = threshold + 1
            stop = True

        for elt in interactions_count:
            if c_tot[threshold] != 0:
                interactions_count[elt][0][threshold] = interactions_count[elt][0][threshold]/c_tot[threshold]*100
    lines = []
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.set_title('Cutoff '+str(cutoff)+' Å')
    ax.set_xlabel('Threshold')
    ax.set_ylabel(ylabel='Percentage of interaction highlighted by the perturbation network')
    for elt in interactions_count:
        color = next(colors)
        lines.append(mlines.Line2D([], [], color=color, label=interactions_count[elt][1]))
        if cutoff == 3:
            ax.plot(np.arange(0, 1.1, 0.1).tolist(), interactions_count[elt][0], c=color)
        elif stop:
            ax.plot(L_cutoffs[:stopping_point], interactions_count[elt][0][:stopping_point], c=color)
        else:
            ax.plot(L_cutoffs, interactions_count[elt][0], c=color)

    f.legend(handles=lines, loc='upper left', fontsize=8)
    plt.savefig('/home/agheerae/results/type_of_interactions_'+str(cutoff)+'.pdf')








    

