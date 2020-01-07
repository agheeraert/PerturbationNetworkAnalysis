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
from Bio.PDB.Polypeptide import aa3
from math import sqrt, ceil

# NW1 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_apo.p", "rb"))
# NW2 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_prfar.p", "rb"))

# L_aa = set([aa for elt in NW1 for aa in elt])
# NW1_c, NW2_c = [0]*len(L_aa), [0]*len(L_aa)
# for aa_pos, node in enumerate(L_aa):
#     if node[-1] == "H":
#         node_str = node[0]+str(int(node[1:-2])-253)+node[-1]
#     else:
#         node_str = node
#     for i in range(len(NW1)):
#         if node in NW1[i]:
#             NW1_c[aa_pos]+=1
#         if node in NW2[i]:
#             NW2_c[aa_pos]+=1

# D_NW = [a-b for a, b in zip(NW2_c, NW1_c)]
# f = plt.figure()
# plt.subplot(121)
# plt.bar(range(len(NW1_c)), NW1_c, color = 'b', width=1)
# plt.subplot(122)
# plt.bar(range(len(NW2_c)), NW2_c, color = 'r', width=1)
# plt.savefig('/home/aghee/results/long_range_evolution/count_nodes_in_LRN.png')
# plt.close()

# f = plt.figure()
# plt.subplot(121)
# plt.bar(range(len(NW1_c)), sorted(NW1_c, reverse=True), color = 'b', width=1)
# plt.subplot(122)
# plt.bar(range(len(NW2_c)), sorted(NW2_c, reverse=True), color = 'r', width=1)
# plt.savefig('/home/aghee/results/long_range_evolution/ranked_nodes_in_LRN.png')
# plt.close()

# f = plt.figure
# plt.bar(range(len(D_NW)), sorted(D_NW, reverse=True), color = 'black', width=1)
# plt.savefig('/home/aghee/results/long_range_evolution/delta_nodes_in_LRN.png')
# plt.close()

# L_aa = list(L_aa)
# for i, node in enumerate(L_aa):
#     if node[-1] == 'H':
#         L_aa[i] = node[0]+str(int(node[1:-2])-253)+":"+node[-1]

# print(sorted(zip(D_NW, L_aa)))

NW1_c, NW2_c = [0]*454, [0]*454
for i in range(1,5):
    NW1 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_apo_"+str(i)+".p", "rb"))
    NW2 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_prfar_"+str(i)+".p", "rb"))
    L_aa = set([aa for elt in NW1 for aa in elt])
    for aa_pos, node in enumerate(L_aa):
        if node[-1] == "H":
            node_str = node[0]+str(int(node[1:-2])-253)+node[-1]
        else:
            node_str = node
        for i in range(len(NW1)):
            if node in NW1[i]:
                NW1_c[aa_pos]+=1
            if node in NW2[i]:
                NW2_c[aa_pos]+=1

D_NW = [a-b for a, b in zip(NW2_c, NW1_c)]
f = plt.figure()
plt.subplot(121)
plt.bar(range(len(NW1_c)), NW1_c, color = 'b', width=1)
plt.subplot(122)
plt.bar(range(len(NW2_c)), NW2_c, color = 'r', width=1)
plt.savefig('/home/aghee/results/long_range_evolution/count_nodes_in_LRN_global.png')
plt.close()

f = plt.figure()
plt.subplot(121)
plt.bar(range(len(NW1_c)), sorted(NW1_c, reverse=True), color = 'b', width=1)
plt.subplot(122)
plt.bar(range(len(NW2_c)), sorted(NW2_c, reverse=True), color = 'r', width=1)
plt.savefig('/home/aghee/results/long_range_evolution/ranked_nodes_in_LRN_global.png')
plt.close()

f = plt.figure
plt.bar(range(len(D_NW)), sorted(D_NW, reverse=True), color = 'black', width=1)
plt.savefig('/home/aghee/results/long_range_evolution/delta_nodes_in_LRN_global.png')
plt.close()

L_aa = list(L_aa)
for i, node in enumerate(L_aa):
    if node[-1] == 'H':
        L_aa[i] = node[0]+str(int(node[1:-2])-253)+":"+node[-1]

print(sorted(zip(D_NW, L_aa))[:11], sorted(zip(D_NW, L_aa))[-10:])
print(sorted(zip(NW1_c, L_aa))[:11])
print(sorted(zip(NW2_c, L_aa))[:11])