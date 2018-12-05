from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

path = 'data/sim1/apo/1frame_1.pdb'

structure = PDBParser().get_structure('X', path)[0]
three2one = dict(zip(aa3, aa1))
hisH, hisF = {}, {}

for residue in structure.get_residues():
    if residue.parent.id == 'H':
        hisH[residue.id[1]] = three2one[residue.resname]
    if residue.parent.id == 'F':
        hisF[residue.id[1]] = three2one[residue.resname]

loop1 = {"chain": 'F', "start": 17, "end": 30}
b2 = {"chain": 'F', "start": 46, "end": 51}
a2 = {"chain": 'F', "start": 56, "end": 72}
a3 = {"chain": 'F', "start": 85, "end": 95}
oloop_and_a1 = {"chain": 'H', "start": 8, "end": 26}
pgvg = {"chain": 'H', "start": 49, "end": 52}
c84 = {"chain": 'H', "start": 84, "end": 84}
allopathway = [loop1, b2, a2, a3, oloop_and_a1, pgvg, c84]

allopathway_nodes = []
for elt in allopathway:
    for i in range(elt["start"], elt["end"]+1):
        allopathway_nodes.append(hisF[i]+str(i)+':'+elt["chain"])
b1F = {"chain": 'F', "start": 1, "end": 16}
a1F = {"chain": 'F', "start": 30, "end": 45}
b3F = {"chain": 'F', "start": 73, "end": 84}
b4F = {"chain": 'F', "start": 95, "end": 102}


if __name__ == '__main__':
    print(allopathway_nodes, len(allopathway_nodes)/500)

