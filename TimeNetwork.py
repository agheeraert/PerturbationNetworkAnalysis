from CreatePerturbationNetwork import PerturbationNetwork
from DrawNetwork import DrawNetwork
from multiprocessing import Pool
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import itertools
import os
from os.path import join, dirname, isdir, basename
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

num = itertools.cycle((1, 2))

class TimePertNet(PerturbationNetwork):
    """Creates the Time Perturbation Network"""
    def __init__(self, path1, path2, out_dir, cutoff=5):
        super().__init__(cutoff=cutoff)