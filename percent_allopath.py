import numpy as np
from allopathway import allopathway_nodes
import os
import networkx as nx
import pandas as pd

base_folder = '/home/agheerae/results/cutoff/'
output1 = '/home/agheerae/results/cutoff/allosteric_path1.xlsx'
output2 = '/home/agheerae/results/cutoff/allosteric_path2.xlsx'
output3 = '/home/agheerae/results/cutoff/allosteric_path_ref.xlsx'

res1 = np.zeros([7, 34])
res2 = np.zeros([7, 34])
res1[:,0] = np.asarray(list(range(3, 10)))
res2[:,0] = np.asarray(list(range(3, 10)))
superbonus1, superbonus2, superbonust = {}, {}, {}
for _dir in os.listdir(base_folder):
    d = os.path.join(base_folder, _dir)
    if os.path.isdir(d):
        for _file in os.listdir(d):
            if _file[-2:] == '.p':
                cutoff = int(_file[0])
                if _file[3] == '-':
                    thresh = float(_file[2:5].replace('-', '.'))
                else:
                    try:
                            thresh = int(_file[2:4])
                    except ValueError:
                            thresh = int(_file[2])
                net = nx.read_gpickle(os.path.join(d, _file))
                c1 = 0
                for edge in net.edges():
                    if edge[0][1:] in allopathway_nodes and edge[1][1:] in allopathway_nodes:
                        c1 +=1
                if isinstance(thresh, int):
                    res1[cutoff-3, thresh] = float(c1)/float(len(net.edges()))
                    res2[cutoff-3, thresh] = len(net.edges())
                else:
                    if cutoff in superbonus1:   
                        superbonus1[cutoff].append(c1/len(net.edges()))
                        superbonust[cutoff].append(thresh)
                        superbonus2[cutoff].append(len(net.edges()))
                    else: 
                        superbonus1[cutoff]=[c1/len(net.edges())]
                        superbonust[cutoff]=[ thresh ]
                        superbonus2[cutoff]=[len(net.edges())]
first = True
for elt in superbonus1:
    _array = np.stack([np.asarray(superbonust[elt]), np.asarray(superbonus1[elt]), np.asarray(superbonus2[elt])])
    if not first:
        array = np.concatenate([array, _array], axis=1)
    else:
        array = _array
        first = False

df1 = pd.DataFrame(res1)
df2 = pd.DataFrame(res2)
df3 = pd.DataFrame(array)
df1.to_excel(output1)
df2.to_excel(output2)
df3.to_excel(output3)





