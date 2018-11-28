import argparse
import shutil
import os.path as o

parser = argparse.ArgumentParser(description='Permutes H and F for my wrong files')
parser.add_argument('path',  type=str,
                     help='file to permute')
args = parser.parse_args()

dire = o.dirname(args.path)
temp = o.join(dire, 'temp.pdb')

tempfile = open(temp, 'w')

with open(args.path, 'r') as fichier:
    for line in fichier:
        if len(line) >= 21:
            if line[23] in ['2','3','4']:
                nb = int(line[23:26]) - 253
                if nb > 0:
                    if len(str(nb)) == 1:
                        tempfile.write(line[:23]+'  '+str(nb)+line[26:])
                    elif len(str(nb)) == 2:
                        tempfile.write(line[:23]+' '+str(nb)+line[26:])
                    else:
                        tempfile.write(line[:23]+str(nb)+line[26:])
                else:
                    tempfile.write(line)
            else:
                tempfile.write(line)
        else:
            tempfile.write(line)
shutil.move(temp, args.path)