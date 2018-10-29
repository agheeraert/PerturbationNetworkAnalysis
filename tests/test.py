def alnParser(id1, id2, alnpath):
    with open(alnpath, 'r') as aln_file:
        chain1, chain2, one2two, two2one = None, None, None, None
        for line in aln_file:
            words = line.split()
            if words:
                if words[0:2] == ['Chain', '1:']:
                    chain1 = [aa + str(i) + ':' + id1 for i, aa in enumerate(list(words[3]), int(words[2]))]
                if words[0:2] == ['Chain', '2:']:
                    chain2 = [aa + str(i) + ':' + id2 for i, aa in enumerate(list(words[3]), int(words[2]))]
            if chain1 and chain2:
                if not (one2two and two2one):
                    one2two = dict(zip(chain1, chain2))
                    two2one = dict(zip(chain2, chain1))
                if one2two and two2one:
                    one2two.update(dict(zip(chain1, chain2)).items())
                    two2one.update(dict(zip(chain2, chain1)).items())
    return one2two, two2one

if __name__ == '__main__':
    print(alnParser(id1='A', id2='B', alnpath='/home/hgheerae/Python/PerturbationNetworkAnalysis/tests/test.aln'))            