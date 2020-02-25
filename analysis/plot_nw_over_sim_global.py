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

NW_apo, NW_prfar = [], []
for i in range(1,5):
    NW_apo.append(pkl.load(open("/home/aghee/results/long_range_evolution/NW_apo_"+str(i)+".p", "rb")))
    NW_prfar.append(pkl.load(open("/home/aghee/results/long_range_evolution/NW_prfar_"+str(i)+".p", "rb")))

L_aa = set([aa for elt in NW_apo[0] for aa in elt])

mean_nw = np.zeros((1000, 4, 2))
for frame in range(0, 1000):
    for sim in range(4):
        _mean_nw_apo, _mean_nw_prfar = 0, 0
        for node in L_aa:
            if node in NW_apo[sim][frame]:
                _mean_nw_apo += NW_apo[sim][frame][node]
            if node in NW_prfar[sim][frame]:
                _mean_nw_prfar += NW_prfar[sim][frame][node]
        mean_nw[frame, sim, 0] = _mean_nw_apo/len(L_aa)
        mean_nw[frame, sim, 1] = _mean_nw_prfar/len(L_aa)

f = plt.figure()
for sim in range(4):
    plt.subplot(2,2,sim+1)
    plt.plot(range(1000), mean_nw[:,sim,0], color = 'b')
    plt.plot(range(1000), mean_nw[:,sim,1], color = 'r')
    plt.title("Simulation "+str(sim+1))
lines = []
lines.append(mlines.Line2D([], [], color='b', label='APO'))
lines.append(mlines.Line2D([], [], color='r', label='PRFAR'))
f.legend(handles=lines, fontsize=8, bbox_to_anchor=(0.65, 1), ncol=2, fancybox=True, shadow=True)
plt.tight_layout()
plt.savefig('/home/aghee/results/long_range_evolution/global.png')
plt.close()
