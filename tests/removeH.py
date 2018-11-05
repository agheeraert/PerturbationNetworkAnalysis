import argparse

parser = argparse.ArgumentParser(description='file paths')
parser.add_argument('filepath',  type=str,
                    help='input file')
parser.add_argument('out',  type=str,
                    help='output file')
args = parser.parse_args()

output = open(args.out, 'w')

with open(args.filepath, 'r') as fichier:
    for line in fichier:
        if len(line.split()) < 2 or line.split()[2][0] != 'H':
            output.write(line)