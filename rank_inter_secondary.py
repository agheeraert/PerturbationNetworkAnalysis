import networkx as nx
import numpy as np
from os.path import dirname
from os.path import join as jn
import matplotlib.pyplot as plt
import itertools
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from math import sqrt, copysign
import warnings
import argparse
import pandas as pd
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
three2one = dict(zip(aa3, aa1))

parser = argparse.ArgumentParser(description='Load a graph and ranks the interactions between secondary structures.')
parser.add_argument('path',  type=str,
                     help='file to load')
parser.add_argument('output',  type=str,
                     help='output path')

args = parser.parse_args()

net = nx.read_gpickle(args.path)

secondary_structure = { (3, 12): "f-b1",
(12, 28): "loop1",
(29, 42): "f-a1",
(45, 51): "f-b2",
(59, 71): "f-a2",
(75, 80): "f-b3",
(82, 95): "f-a3",
(98, 100): "f-b4",
(102, 106): "f-a4p",
(109, 117): "f-a4",
(121, 124): "f-a4pp",
(125, 133): "f-b5",
(135, 140): "f-b5p",
(154, 163): "f-a5",
(166, 169): "f-b6",
(170, 174): "f-a6p",
(182, 185): "f-a6",
(195, 200): "f-b7",
(207, 217): "f-a7",
(218, 223): "f-b8",
(231, 237): "f-a8p",
(239, 244): "f-a8",
(253, 259): "h-b1",
(260, 265): "omega-loop",
(266, 276): "h-a1",
(281, 288): "h-b2",
(295, 300): "h-b3",
(301, 304): "PGVG",
(305, 316): "h-a2",
(318, 329): "h-a2p",
(330, 336): "h-b4",
(337, 341): "h-a3",
(345, 347): "h-b5",
(359, 363): "h-b6",
(371, 382): "h-b7",
(386, 399): "h-b8",
(402, 408): "h-b9",
(415, 420): "h-b10",
(423, 428): "h-b11",
(434, 449): "h-a4"}

sec_all = {}

for liste in secondary_structure:
    for i in range(liste[0], liste[1]+1):
        sec_all[i] = secondary_structure[(liste[0], liste[1])]

new, count, still = True, 0, False

for j in range(454):
    if j not in sec_all:
        if new:
            count+=1
            new = False
            still = True
        sec_all[j]= 'segment-'+str(count)
    else:
        if still:
            still = False
            new = True
node_to_sec = {}

for node in net.nodes():
    if node[-1] == 'H':
        c = 252
    else:
        c = -1
    if c+int(node[1:-2]) in sec_all:
        node_to_sec[node] = sec_all[c+int(node[1:-2])]
    else:
        node_to_sec[node] = 'unknown'

print(node_to_sec)

c_inter_sec = {}
for u, v in net.edges():
    w = net.get_edge_data(u, v)['weight']
    sec_sec = (node_to_sec[u], node_to_sec[v])
    if sec_sec in c_inter_sec:
        c_inter_sec[sec_sec] +=w
    elif sec_sec[::-1] in c_inter_sec:
        c_inter_sec[sec_sec[::-1]] = w
    else:
        c_inter_sec[sec_sec] = w
# print(sorted(c_inter_sec, key=c_inter_sec.get, reverse=True))
df = pd.DataFrame(c_inter_sec, index=[0]).T.sort_values(0, axis=0, ascending=False)
df.to_excel(args.output)
