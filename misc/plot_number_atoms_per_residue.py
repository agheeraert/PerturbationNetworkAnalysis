import matplotlib.pyplot as plt 

dict_size_letter = {'A':6,
'C':7,
'D':9,
'E':10,
'F':12,
'G':5,
'H':11,
'I':9,
'K':10,
'L':9,
'M':9,
'N':9,
'P':8,
'Q':10,
'R':12,
'S':7,
'T':8,
'V':8,
'W':15,
'Y':13 } 

ordered = [ (key , value) for (key, value) in sorted(dict_size_letter.items() ,  key=lambda x: x[1]  ) ]
plt.bar([elt[0] for elt in ordered], [elt[1] for elt in ordered], color='cornflowerblue')
plt.title('Ranking of amino acid in number of heavy atoms')
plt.xlabel('Amino acid')
plt.ylabel('Number of heavy atoms')
plt.savefig('number_of_heavy_atoms_per_residue.png')




