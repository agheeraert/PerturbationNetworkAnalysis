import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from tqdm import tqdm
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import argparse

from CreateNetwork import three2one, AANetwork

parser = argparse.ArgumentParser(description='Get a network and outputs only the interface(s)')
parser.add_argument('f',  type=str, nargs='+',
                    help='List of network files')
parser.add_argument('-v', action='store_true', default=False,
                    help='Outputs table with interface residues, node and ')
args = parser.parse_args()

for f in args.f:
    net = AANetwork()
    net.net = nx.read_gpickle(f)
    if args.v:
        net.get_interface(f.replace('.p', '.xlsx'))
    else:
        net.get_interface()
    assert '.p' in f, print('Should use the .p extension')
    net.save(f.replace('.p', '_i.p'))
        

