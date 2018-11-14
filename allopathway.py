loop1 = {"chain": 'F', "start": 17, "end": 30}
b2 = {"chain": 'F', "start": 46, "end": 51}
a2 = {"chain": 'F', "start": 56, "end": 72}
a3 = {"chain": 'F', "start": 85, "end": 95}
oloop_and_a1 = {"chain": 'H', "start": 8, "end": 26}
pgvg = {"chain": 'H', "start": 49, "end": 72}
c84 = {"chain": 'H', "start": 84, "end": 84}
allopathway = [loop1, b2, a2, a3, oloop_and_a1, pgvg, c84]

allopathway_nodes = []
for elt in allopathway:
    for i in range(elt["start"], elt["end"]+1):
        allopathway_nodes.append(i+':'+elt["chain"])