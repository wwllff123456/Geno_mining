import Align


def test():
    d = open('example_data\\DNA.txt', 'r')
    g = open('example_data\\gRNA.txt', 'r')
    r = open('example_data\\text.txt', 'w')

    dna = d.read().replace('\n', '')

    print(len(dna))

    for i in range(21, len(dna)-2):
        if dna[i:i+2] == 'GG':
            r.write(dna[i-21:i+2])
            r.write('\n')

    dna = dna.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

    for i in range(21, len(dna)-2):
        if dna[i:i+2] == 'GG':
            r.write(dna[i-21:i+2])
            r.write('\n')

    d.close()
    g.close()
    r.close()


if __name__ == '__main__':
    test()