from os import listdir
from os.path import join as jn
import networkx as nx

BASE_FOLDER = '/home/agheerae/results/avg_sims/pertnet/'

dict_residues = {['R', 'H', 'K']: '+',
                ['D', 'E']: '-',
                ['S', 'T', 'N', 'Q']: 'p',
                ['G', 'A', 'P', 'Y']: '?',
                ['V', 'I', 'L', 'M', 'F', 'W', 'C']: 'h'}

for liste in dict_residues:
    for elt in liste:
        dict_residues[elt] = dict_residues[liste]

dict_interactions = {('+', '-'): 'sb',
                    ('+', '+'): 'un_sb',
                    ('-', '-'): 'un_sb',
                    ('h', 'h'): 'h',
                    ('p', 'p'): 'p'
}

interactions_count = {'sb': 0,
                     'un_sb': 0,
                     'h': 0,
                     'p': 0,
                     '?': 0
}

for cutoff in range(3, 10):
    folder = jn(BASE_FOLDER, 'cutoff_'+str(cutoff))
    if cutoff = 3:
        L_cutoffs = range(0, 11)
    else:
        L_cutoffs = range(0, 30)
    for threshold in L_cutoffs:
        if cutoff = 3:
            if threshold in [0, 10]:
                string = str(threshold)[0]
            else:
                string = '0-'+str(threshold)
        else:
            string = str(threshold)
        try:
            net = nx.read_gpickle(jn(folder, str(cutoff)+'_'+string+'.p'))
            for u, v in A.edges():
                residue_residue = (dict_residues[u[0]], dict_residues[v[0]])
                if residue_residue in dict_interactions:
                    interaction = dict_interactions[residue_residue]
                elif residue_residue[::-1] in dict_interactions:
                    interaction = dict_interactions[residue_residue]
                else:
                    interaction = '?'
                interaction_count[interaction]+=1




    

