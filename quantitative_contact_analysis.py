from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from os import listdir
from os.path import join as jn
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

class args:
    base_residue = 'K19:F'
    residues_to_watch = 'H228:F'
    method = 'crd'
    frames_apo_paths = [jn('data/H/apo/', path) for i, path in enumerate(listdir('data/H/apo/')) if (path[0] == '1' and i%10==0)]
    frames_holo_paths = [jn('data/H/prfar/', path) for i, path in enumerate(listdir('data/H/prfar/')) if (path[0] == '1' and i%10==0)]
    output = '/home/agheerae/results/K19:F-H228:f_10_.svg'

class CationRingDistance():
    def __init__(self, cation, ring):
        cation2atoms = {'K': ['HZ1', 'HZ2', 'HZ3']}
        ring2atoms = {'H': ['CG', 'ND1', 'CE1', 'NE2', 'CD2']}
        self.cation_atoms = cation2atoms[cation[0]]
        self.cation_position = int((cation.split(':')[0])[1:]) 
        self.cation_chain = cation.split(':')[1]
        self.ring_atoms = ring2atoms[ring[0]] 
        self.ring_position = int((ring.split(':')[0])[1:])
        self.ring_chain = ring.split(':')[1]

    def distance(self, pdb):
        structure = PDBParser().get_structure('X', pdb)[0]
        L_cation_coords, L_ring_coords = [], []

        for elt in self.cation_atoms:
            L_cation_coords.append(structure[self.cation_chain][(' ', self.cation_position, ' ')][elt].coord)
        cation_coord = np.mean(L_cation_coords, axis=0)        
        for elt in self.ring_atoms:
            L_ring_coords.append(structure[self.ring_chain][(' ', self.ring_position, ' ')][elt].coord)
        ring_coord = np.mean(L_ring_coords, axis=0)
        return norm(ring_coord-cation_coord)

method2class = {'crd': CationRingDistance}

def traj_distance(method, base, watch, frames):
    Distance = method2class[method](base, watch)
    distance_array = np.zeros(len(frames))
    for i, filepath in enumerate(frames):
        distance_array[i] = Distance.distance(filepath)
    return distance_array

f = plt.figure()
plt.title('Evolution of the distance between '+ args.base_residue+  'and ' + args.residues_to_watch + ' during MD simulation')
plt.xlabel('Time (in ns)')
plt.ylabel('Distance between ' + args.base_residue + 'and '+ args.residues_to_watch +' (in Å)')
X = np.arange(0, (len(args.frames_apo_paths)/10), 0.1)
Y1 = traj_distance(args.method, args.base_residue, args.residues_to_watch, args.frames_apo_paths)
Y2 = traj_distance(args.method, args.base_residue, args.residues_to_watch, args.frames_holo_paths)
plt.plot(X, Y1, c='g')
plt.plot(X, Y2, c='r')
plt.savefig(args.output)


	
	
