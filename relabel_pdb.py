import argparse
import shutil
from os.path import join as jn, dirname
from os import remove as rm
parser = argparse.ArgumentParser(description='Relabel correctly PDB frames')
parser.add_argument('f',  type=str, nargs='+',
                     help='files to relabel')
args = parser.parse_args()

offset = 1
chain_separation = 253

for frame in args.f:
    wd = dirname(frame)
    temp = jn(wd, 'temp.pdb')
    tempfile = open(temp, 'w')
    with open(frame, 'r') as fichier:
        for line in fichier:
            if line.split()[0] == 'ATOM':
                aa_id = int(line.split()[5])
                if aa_id < chain_separation:
                    chain = 'F'
                    new_id = str(aa_id + 1)
                else:
                    chain = 'H'
                    new_id = str(aa_id - chain_separation + 1)
                if len(new_id) == len(str(aa_id)):
                    new_id = ' %s '%new_id
                elif len(new_id) > len(str(aa_id)):
                    new_id = new_id+' '
                else:
                    n = len(str(aa_id))-len(new_id)+1
                    new_id = ' '*n+new_id+' '
                line = line.replace(' %s  ' %str(aa_id), new_id).replace(' X ', ' %s ' %chain)
            tempfile.write(line)
    shutil.move(temp, frame)
