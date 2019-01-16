from os import listdir
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import pickle as pkl
import pandas as pd

BASE_FOLDER = '/home/agheerae/results/sidechain_H/pertnet/'
OUT_FOLDER = '/home/agheerae/results/type_of_interactions/sidechain_H/'
colors = itertools.cycle(('r', 'g', 'b', 'olive', 'purple', 'navy'))   # 'purple',  between olive and navy
tuple_residues = [(['R', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y', 'H'], 'p'),
                (['G', 'A', 'P'], '?'),
                (['V', 'I', 'L', 'M', 'F', 'W', 'C'], 'h')]
aa_list = [aa for tup in tuple_residues for aa in tup[0]]
q = len(aa_list)
aa_to_id = dict(zip(aa_list, range(q)))
n_top_unknown = 10

dict_residues = {}
for tup in tuple_residues:
    for elt in tup[0]:
        dict_residues[elt] = tup[1]

def update_dico(dico, tup):
    if tup in dico:
        dico[tup]+=1
    elif tup[::-1] in dico:
        dico[tup[::-1]]+=1
    else:
        dico[tup] = 1
    return dico

dict_interactions = {('+', '-'): 'sb',
                    ('+', '+'): 'un_sb',
                    ('-', '-'): 'un_sb',
                    ('h', 'h'): 'h',
                    ('p', 'p'): 'p'
}

for cutoff in range(5, 6):
    folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
    last_cutoff = int((len(listdir(folder))-1)/2)
    L_threshs = range(0, last_cutoff+1)
    unknown_pairs = []
    known_pairs = []
    interactions_count = {'sb': ((L_threshs[-1]+1)*[0], "Salt bridge"),
        'un_sb': ((L_threshs[-1]+1)*[0], "Unknown salt bridge"),
        'h': ((L_threshs[-1]+1)*[0], "Hydrophobic interaction"),
        'p': ((L_threshs[-1]+1)*[0], "Polar-polar interaction"),
        '?': ((L_threshs[-1]+1)*[0], "Unknown interaction"),
        'c': ((L_threshs[-1]+1)*[0], "Residues covalently bound"),
                }
    c_tot = (L_threshs[-1]+1)*[0]
    for threshold in L_threshs:
        _unknown_pairs = {}
        _known_pairs = {}
        if str(cutoff)+'_'+str(threshold)+'.p' in listdir(folder):
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+str(threshold)+'.p'))
            for u, v in net.edges():
                if abs(int(u[1:].split(':')[0]) - int(v[1:].split(':')[0])) == 1:
                    interaction = 'c'
                else:
                    residue_residue = (dict_residues[u[0]], dict_residues[v[0]])
                    if residue_residue in dict_interactions:
                        interaction = dict_interactions[residue_residue]
                        _known_pairs = update_dico(_known_pairs, (u[0], v[0]))
                    elif residue_residue[::-1] in dict_interactions:
                        interaction = dict_interactions[residue_residue[::-1]]
                        _known_pairs = update_dico(_known_pairs, (u[0], v[0]))
                    else:
                        interaction = '?'
                        _unknown_pairs = update_dico(_unknown_pairs, (u[0], v[0]))
                interactions_count[interaction][0][threshold]+=1
                c_tot[threshold]+=1
            unknown_pairs.append(_unknown_pairs)
            known_pairs.append(_known_pairs)
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
        ax.plot(L_threshs, interactions_count[elt][0], c=color)

    f.legend(handles=lines, loc='upper left', fontsize=8)
    plt.savefig(jn(OUT_FOLDER, 'type_of_interactions_'+str(cutoff)+'.svg'))
    
    L_top_unknown = []
    for i, elt in enumerate(unknown_pairs):
        f = plt.figure()
        mat = np.zeros([q, q])
        for u, v in elt:
            mat[(aa_to_id[u], aa_to_id[v])] = elt[(u, v)]
            mat[(aa_to_id[v], aa_to_id[u])] = elt[(u, v)]
        plt.imshow(mat, cmap = 'Reds')
        plt.colorbar()
        plt.xticks(range(q), aa_list)
        plt.yticks(range(q), aa_list)
        plt.savefig(jn(OUT_FOLDER, 'unknown_pairs_'+str(cutoff)+'_'+str(i)+'.svg'))
        f = plt.figure()
        for u, v in known_pairs[i]:
            mat[(aa_to_id[u], aa_to_id[v])] += known_pairs[i][(u, v)]
            mat[(aa_to_id[v], aa_to_id[u])] += known_pairs[i][(u, v)]
        plt.imshow(mat, cmap = 'Reds')
        plt.colorbar()
        plt.xticks(range(q), aa_list)
        plt.yticks(range(q), aa_list)
        plt.savefig(jn(OUT_FOLDER, 'all_pairs_'+str(cutoff)+'_'+str(i)+'.svg'))
        _L_top_unknown = []
        for j, tup in enumerate(sorted(unknown_pairs[i], key=unknown_pairs[i].get, reverse=True)[:n_top_unknown]):
            _L_top_unknown.extend([tup, unknown_pairs[i][tup]])
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+str(i)+'.p'))
            avg = 0
            for u, v in net.edges():
                if (u[0], v[0]) == tup or (v[0], u[0]) == tup:
                    avg += net.get_edge_data(u, v)['weight']
                    print(u, v, net.get_edge_data(u, v)['weight'])
            avg /= unknown_pairs[i][tup]
            _L_top_unknown.append(avg)
        L_top_unknown.append(_L_top_unknown)
    df = pd.DataFrame(L_top_unknown, index=range(last_cutoff+1))
    df.to_csv(OUT_FOLDER+'top_unknown.csv')