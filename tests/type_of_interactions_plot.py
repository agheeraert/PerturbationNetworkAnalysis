from os import listdir
from os.path import join as jn
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np
import itertools
import matplotlib.lines as mlines
import matplotlib
import pickle as pkl
import pandas as pd

BASE_FOLDER = '/home/agheerae/results/avg_sims/pertnet/'
OUT_FOLDER = '/home/agheerae/results/type_of_interactions/without_H/'
colors = itertools.cycle(('r', 'g', 'b', 'olive', 'purple', 'navy'))   # 'purple',  between olive and navy
tuple_residues = [(['R', 'K'], '+'),
                (['D', 'E'], '-'),
                (['S', 'T', 'N', 'Q', 'Y', 'H'], 'p'),
                (['I', 'L', 'V', 'M', 'F', 'W', 'C', 'P', 'G', 'A'], 'h')]
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
# matplotlib.rcParams['font.family'] = "serif"

for cutoff in range(5, 6):
    folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
    last_cutoff = int((len(listdir(folder))-1)/3)
    L_threshs = range(0, 9)
    unknown_pairs = []
    known_pairs = []
    interactions_count = {'sb': ((L_threshs[-1]+1)*[0], "Salt bridge"),
        'un_sb': ((L_threshs[-1]+1)*[0], "Same charge"),
        'h': ((L_threshs[-1]+1)*[0], "Hydrophobic"),
        'p': ((L_threshs[-1]+1)*[0], "Polar"),
        '?': ((L_threshs[-1]+1)*[0], "Non-specific"),
        'c': ((L_threshs[-1]+1)*[0], "Covalently bound"),
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
    ax.set_xlabel('Threshold')
    ax.set_ylabel(ylabel='Percentage of interactions in the perturbation network')
    for elt in interactions_count:
        color = next(colors)
        lines.append(mlines.Line2D([], [], color=color, label=interactions_count[elt][1]))
        ax.plot(L_threshs, interactions_count[elt][0], c=color)
    ax2 = ax.twinx()
    ax2.semilogy(L_threshs, c_tot, color='black')
    lines.append(mlines.Line2D([], [], color='black', label='Total'))
    ax2.set_ylabel(ylabel='Number of interactions in the perturbation network')
    ax.grid(linestyle='--')
    f.legend(handles=lines, loc='upper center', fontsize=8, bbox_to_anchor=(0.4, 1),
          ncol=2, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(jn(OUT_FOLDER, 'type_of_interactions_article2_'+str(cutoff)+'.png'))
    
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
        plt.savefig(jn(OUT_FOLDER, 'unknown_pairs_'+str(cutoff)+'_'+str(i)+'.png'))
        f = plt.figure()
        for u, v in known_pairs[i]:
            mat[(aa_to_id[u], aa_to_id[v])] += known_pairs[i][(u, v)]
            mat[(aa_to_id[v], aa_to_id[u])] += known_pairs[i][(u, v)]
        plt.imshow(mat, cmap = 'Reds')
        plt.colorbar()
        plt.xticks(range(q), aa_list, family='serif')
        plt.yticks(range(q), aa_list, family='serif')
        plt.savefig(jn(OUT_FOLDER, 'all_pairs_'+str(cutoff)+'_'+str(i)+'.png'))
        _L_top_unknown = []
        for j, tup in enumerate(sorted(unknown_pairs[i], key=unknown_pairs[i].get, reverse=True)[:n_top_unknown]):
            _L_top_unknown.extend([tup, unknown_pairs[i][tup]])
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+str(i)+'.p'))
            avg = 0
            for u, v in net.edges():
                if (u[0], v[0]) == tup or (v[0], u[0]) == tup:
                    avg += net.get_edge_data(u, v)['weight']
            avg /= unknown_pairs[i][tup]
            _L_top_unknown.append(avg)
        L_top_unknown.append(_L_top_unknown)
    df = pd.DataFrame(L_top_unknown, index=range(last_cutoff+1))
    df.to_csv(OUT_FOLDER+'top_unknown_article_2.csv')
