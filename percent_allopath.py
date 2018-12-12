import numpy as np
from allopathway import allopathway_nodes
import os
import networkx as nx
import pandas as pd

base_folder = '/home/agheerae/results/timewin_500/pertnet/'
out_edge1 = os.path.join(base_folder, 'allosteric_path1_edges.xlsx')
out_edge2 = os.path.join(base_folder, 'allosteric_path1s_edges.xlsx')
out_edge3 = os.path.join(base_folder, 'allosteric_path2_edges.xlsx')
out_node1 = os.path.join(base_folder, 'allosteric_path1_nodes.xlsx')
out_node2 = os.path.join(base_folder, 'allosteric_path1s_nodes.xlsx')
out_node3 = os.path.join(base_folder, 'allosteric_path2_nodes.xlsx')

res1 = np.zeros([7, 100])
res2 = np.zeros([7, 100])
res3 = np.zeros([7, 100])
res4 = np.zeros([7, 100])
res1[:,0] = np.asarray(list(range(3, 10)))
res2[:,0] = np.asarray(list(range(3, 10)))
res3[:,0] = np.asarray(list(range(3, 10)))
res4[:,0] = np.asarray(list(range(3, 10)))

superbonus1, superbonus2, superbonust, superbonus3, superbonus4, superbonustt = {}, {}, {}, {}, {}, {}
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
                c1, c2 = 0, 0
                for edge in net.edges():
                    if edge[0] in allopathway_nodes and edge[1] in allopathway_nodes:
                        c1 +=1
                for node in net.nodes():
                    if node in allopathway_nodes:
                        c2 +=1
                if isinstance(thresh, int):
                    res1[cutoff-3, thresh] = float(c1)/float(len(net.edges()))
                    res2[cutoff-3, thresh] = len(net.edges())
                    res3[cutoff-3, thresh] = float(c2)/float(len(net.nodes()))
                    res4[cutoff-3, thresh] = len(net.nodes())
                else:
                    if cutoff in superbonus1:   
                        superbonus1[cutoff].append(c1/len(net.edges()))
                        superbonust[cutoff].append(thresh)
                        superbonus2[cutoff].append(len(net.edges()))
                        superbonus3[cutoff].append(c2/len(net.nodes()))
                        superbonustt[cutoff].append(thresh)
                        superbonus4[cutoff].append(len(net.nodes()))
                    else: 
                        superbonus1[cutoff]=[c1/len(net.edges())]
                        superbonust[cutoff]=[ thresh ]
                        superbonus2[cutoff]=[len(net.edges())]
                        superbonus3[cutoff]=[c2/len(net.nodes())]
                        superbonustt[cutoff]=[ thresh ]
                        superbonus4[cutoff]=[len(net.nodes())]                        
first = True
for elt in superbonus1:
    _array = np.stack([np.asarray(superbonust[elt]), np.asarray(superbonus1[elt]), np.asarray(superbonus2[elt])])
    if not first:
        array = np.concatenate([array, _array], axis=1)
    else:
        array = _array
        first = False

first = True
for elt in superbonus3:
    _array = np.stack([np.asarray(superbonustt[elt]), np.asarray(superbonus3[elt]), np.asarray(superbonus4[elt])])
    if not first:
        array2 = np.concatenate([array, _array], axis=1)
    else:
        array2 = _array
        first = False

df1 = pd.DataFrame(res1)
df2 = pd.DataFrame(res2)
df3 = pd.DataFrame(res3)
df4 = pd.DataFrame(res4)
# df5 = pd.DataFrame(array)
# df6 = pd.DataFrame(array2)

df1.to_excel(out_edge1)
df2.to_excel(out_edge2)
df3.to_excel(out_node1)
df4.to_excel(out_node2)
# df5.to_excel(out_edge3)
# df6.to_excel(out_node3)





