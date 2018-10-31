
filepath = '/home/aria/PerturbationNetworkAnalysis/data/frame1000/frame_1000_prfar.pdb'
output = open('/home/aria/PerturbationNetworkAnalysis/data/frame1000/frame_1000_prfar_noH.pdb', 'w')
with open(filepath, 'r') as fichier:
    for line in fichier:
        if len(line.split()) < 2:
            output.write(line) 
            
        elif line.split()[2][0] != 'H':
            output.write(line)