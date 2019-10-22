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

NW1 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_apo.p", "rb"))
NW2 = pkl.load(open("/home/aghee/results/long_range_evolution/NW_prfar.p", "rb"))
L_aa = set([aa for elt in NW1 for aa in elt])
for node in L_aa:
    if node[-1] == "H":
        node_str = node[0]+str(int(node[1:-2])-253)+node[-1]
    else:
        node_str = node
    f = plt.figure()
    NW1_node, NW2_node = [], []
    for i in range(len(NW1)):
        if node in NW1[i]:
            NW1_node.append(NW1[i][node])
        else:
            NW1_node.append(0)
        if node in NW2[i]:
            NW2_node.append(NW2[i][node])
        else:
            NW2_node.append(0)
    plt.plot(list(range(len(NW1))), NW1_node, color = 'b')
    plt.plot(list(range(len(NW2))), NW2_node, color = 'r')
    lines = []
    lines.append(mlines.Line2D([], [], color='b', label='APO'))
    lines.append(mlines.Line2D([], [], color='r', label='PRFAR'))
    f.legend(handles=lines, fontsize=8, bbox_to_anchor=(0.985, 1), ncol=2, fancybox=True, shadow=True)
    plt.title(r"$\bf{" + node_str + "}$")
    plt.tight_layout()
    plt.savefig('/home/aghee/results/long_range_evolution/'+node_str+'.png')
    plt.close()

L_bottom = ["G202:F", "T142:F", "R133:F", "A224:F", "F227:F", "S225:F", "K19:F", "D11:F", "L50:F", "G20:F", "H228:F"]
L_sb = ["E67:F", "E71:F", "R18:H", "E91:F", "R95:F", "Y136:H", "V248:F", "R249:F", "L250:F", "R22:H", "Q72:F", "R187:H", "D74:F", "N247:F", "W123:H", "I73:F", "I75:F", "K184:H", "M14:H"]
L_up = ["N12:H", "N15:H", "P49:H", "G50:H", "V51:H", "G52:H", "G9:H", "P10:H", "G11:H"]
L_unfolding = ["H53:H", "F54:H", "E56:H", "G57:H", "R59:H", "R60:H", "G9:H", "G57:H"]
L_interface = ["M121:H", "R5:F", "K99:F", "E167:F", "D98:F", "K181:H"]


def multiplot(liste, string):
    n = len(liste)
    size = ceil(sqrt(n))
    f = plt.figure()
    for j, node in enumerate(liste):
        plt.subplot(size, size, j+1)
        NW1_node, NW2_node = [], []
        if node[-1] == "H":
            node_str = node[0]+str(int(node[1:-2])+253)+node[-2:]
        else:
            node_str = node
        for i in range(len(NW1)):
            if node_str in NW1[i]:
                NW1_node.append(NW1[i][node_str])
            else:
                NW1_node.append(0)
            if node_str in NW2[i]:
                NW2_node.append(NW2[i][node_str])
            else:
                NW2_node.append(0)
        plt.plot(list(range(len(NW1))), NW1_node, color = 'b')
        plt.plot(list(range(len(NW2))), NW2_node, color = 'r')
        plt.title(r"$\bf{" + node + "}$")
    lines = []
    lines.append(mlines.Line2D([], [], color='b', label='APO'))
    lines.append(mlines.Line2D([], [], color='r', label='PRFAR'))
    f.legend(handles=lines, fontsize=8, loc="lower right", ncol=2, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig('/home/aghee/results/long_range_evolution/'+string+'.svg')
    plt.close()

multiplot(L_bottom, "bottom")
multiplot(L_interface, "interface")
multiplot(L_sb, "sb")
multiplot(L_unfolding, "unfolding")
multiplot(L_up, "up")
