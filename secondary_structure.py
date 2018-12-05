hisF = 'XXXBEEEEEEEEETTEEXXSSSSSXXXXBSXHHHHHHHHHHHTXXEEEEEEXXXSSSHHHHHHHHHHHHHTTXXSXEEEESSXXSHHHHHHHHHTTXSEEEESHHHHTXTHHHHHHHHHHXGGGEEEEEEEEEETTEEEEEETTTTEEEEEEHHHHHHHHHHTTXSEEEEEETTTTTXSSXXXHHHHHHHGGGXXSXEEEESXXXSHHHHHHHHHTTXSEEEESHHHHTTXSXHHHHHHHHHHTTXXBXXXXX'
hisH = 'XEEEEEXXSSSXXHHHHHHHHHHHTTBTTXEEEEESSXXSSXXSEEEEXXXSXSHHHHHHHHHTTXHHHHHHHHHTTXEEEEETHHHHTTSSBBSSSTTXBXXXXSSXEEEEXXXSSXSEEEEEEEEESSSSXXEEEEEEESEEEEXXGGGEEEEEEETTEEEEEEEEETTEEEESSXGGGSHHHHHHHHHHHHHHTTSXX'

secondary = {
    'X': ([],[]), #disorganized
    'B': ([],[]), #bridge
    'S': ([],[]), #bend
    'T': ([],[]), #turn
    'E': ([],[]), #strand
    'G': ([],[]), #3/10 helix
    'H': ([],[]) #helix
}

previous = None
_ = []
for i, structure in enumerate(hisF):
    if structure != previous and len(_) > 0:
        secondary[previous][0].append(_)
        _ = [i+1]
    else:
        _.append(i+1)
    previous = structure
secondary[previous][0].append(_)
previous = None
_ = []

for i, structure in enumerate(hisH):
    if structure != previous and len(_) > 0:
        secondary[previous][1].append(_)
        _ = [i+1]
    else:
        _.append(i+1)
    previous = structure
secondary[previous][1].append(_)
id_to_struct = {}

for structure in secondary:
    for i, elt in enumerate(secondary[structure][0]):
        for j in elt:
            id_to_struct[str(j)+':F'] = structure+str(i+1)

for structure in secondary:
    for i, elt in enumerate(secondary[structure][1]):
        for j in elt:
            id_to_struct[str(j)+':H'] = structure+str(i+1)

print(id_to_struct)