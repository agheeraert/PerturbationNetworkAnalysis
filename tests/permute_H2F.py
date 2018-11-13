import argparse
import shutil
import os.path as o

parser = argparse.ArgumentParser(description='Permutes H and F in all my wrong files...')
parser.add_argument('path',  type=str,
                     help='file to permute')
args = parser.parse_args()

dire = o.dirname(args.path)
temp = o.join(dire, 'temp.pdb')

tempfile = open(temp, 'w')

with open(args.path, 'r') as fichier:
    for line in fichier:
        if len(line) >= 21:
            if line[21] == 'H':
                line = line[:21]+'F'+line[22:]
                tempfile.write(line)
            elif line[21] == 'F':
                line = line[:21]+'H'+line[22:]
                tempfile.write(line)
            else:
                tempfile.write(line)
        else:
            tempfile.write(line)
shutil.move(temp, args.path)