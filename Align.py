import numpy as np
import re
import logging
import datetime

logging.basicConfig(filename='log.log', level=logging.INFO, filemode='w', format='%(message)s')


class DNA:

    def __init__(self):
        self.seq = 'N/A'
        self.name = 'N/A'
        self.file_path = 'N/A'
        logging.info('{}:\n\nStart a new DNA sequence\n'.format(datetime.datetime.now()))

    def set_seq(self, seq: str):
        self.seq = seq
        if len(seq) > 50:
            logging.info('Seq set as {}...\n'.format(seq[:50]))
        else:
            logging.info('Seq set as {}\n'.format(seq))

    def set_name(self, name: str):
        self.name = name
        logging.info('Name set as {}\n'.format(name))

    def set_file_path(self, path: str):
        self.file_path = path
        logging.info('Path set as {}\n'.format(path))

    def __len__(self):
        return len(self.seq)

    def rc(self):
        self.seq = self.seq.upper()
        self.seq = self.seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
        self.seq = self.seq.upper()[::-1]
        logging.info('RC Done.\n')

    def __str__(self):
        return '\nDNA sequence {self.name} at {self.file_path}:\n{self.seq}\n'.format(self=self)


class Guide(DNA):

    pass


class Chromosome(DNA):

    pass


def validate(dna_0: str, dna_1: str, *lists: list, length=20) -> bool:

    valid = True

    for i in lists:
        if len(i) != length:
            valid = False

    if len(dna_0) != length or len(dna_1) != length:
        valid = False

    if not valid:

        print('Not the same length in the two dna or the score lists!!')
        logging.info('Length Error in the two dna or the score lists: ')
        logging.info('\n[{}] {}\n[{}] {}'.format(len(dna_0), dna_0, len(dna_1), dna_1, length))
        for i in lists:
            logging.info('[{}] {}'.format(len(i), i))
        logging.info('Length needs to be {}'.format(length))
        logging.info('Exiting...')
        raise ValueError('Not the same length in the two dna or the score list!!')

    return valid


def linear_compare(dna_0: str, dna_1: str, score_list: list, length=20) -> int:
    """Calculate the score_of_matching for two DNA sequences in a linear way"""

    validate(dna_0, dna_1, score_list, length=20)

    score = 0
    num = 0

    for i in range(len(dna_0)):
        if dna_0[i] == dna_1[i]:
            score += score_list[i]
            num += 1

    logging.info('Compared: \n{}\n{}\n'.format(dna_0, dna_1))
    logging.info('matches: {}/{}\tscore: {}/{}\n'
                 .format(num, len(dna_0), score, sum(score_list)))

    return score


def needleman(dna_0: str, dna_1: str, match_list: list, penalty_list: list, gap_list: list) -> int:
    """Calculate the score_of_matching for two DNA sequences with Needleman-Wunsch algorithm"""

    validate(dna_0, dna_1, match_list, penalty_list, gap_list, length=20)

    penalty = -4
    mis_match = -3
    match = 8

    grid = [[0 for _ in range(21)] for _ in range(21)]
    move = [['' for _ in range(20)] for _ in range(20)]

    grid[0][1] = penalty
    grid[1][0] = penalty
    for i in range(2, 21):
        grid[0][i] = grid[0][i-1] + penalty
        grid[i][0] = grid[i-1][0] + penalty

    for i in range(20):
        for j in range(20):
            '''
            he #i and #j in dna will be i+1 and j+1 in grid because of the adding of first col/row
            this means the diagonal will be [i+1][j+1](now) back to [i][j](up and left), etc.
               
            dna_0 on the left of the grid, dna_1 on top
            '''

            m = mis_match
            if dna_0[i] == dna_1[j]:
                m = match

            # diagonal, vertical, horizontal:
            d, v, h = grid[i][j]+m, grid[i+1][j]+penalty, grid[i][j+1]+penalty

            grid[i+1][j+1] = max(d, v, h)

            if grid[i+1][j+1] == d:
                move[i][j] += 'd'
            if grid[i+1][j+1] == h:
                move[i][j] += 'h'
            if grid[i+1][j+1] == v:
                move[i][j] += 'v'

    way_of_alignment = ['d'*20]

    def align_two(d0, d1, m, d00='', d11='', mm=''):
        if len(d0) == 0 and len(d1) == 0:
            print(d00, d11, mm, sep='\n', end='\n\n')
        try:
            for i in m[-1][-1]:
                if i == 'd':
                    d00 = d0[-1] + d00
                    d11 = d1[-1] + d11
                    d0 = d0[:-1]
                    d1 = d1[:-1]
                    mm += 'd'
                    m = [m[j][:-1] for j in range(len(m)-1)]
                    align_two(d0, d1, m, d00, d11, mm)

                elif i == 'h':
                    d00 = d0[-1] + d00
                    d11 = '-' + d11
                    d0 = d0[:-1]
                    mm += 'h'
                    m = [m[j] for j in range(len(m)-1)]
                    align_two(d0, d1, m, d00, d11, mm)

                elif i == 'v':
                    d00 = '-' + d00
                    d11 = d1[-1] + d11
                    d1 = d1[:-1]
                    mm += 'v'
                    m = [m[j][:-1] for j in range(len(m))]
                    align_two(d0, d1, m, d00, d11, mm)
        except IndexError:
            pass

    def align_two_1(m, mm=''):
        nonlocal way_of_alignment
        if len(m) == 1 or len(m[0]) == 1:
            # print(len(mm), len(m), len(m[0]), mm)
            way_of_alignment.append(mm)
        try:

            if 'd' in m[-1][-1]:
                mm += 'd'
                m = [m[j][:-1] for j in range(len(m) - 1)]
                align_two_1(m, mm)

            else:
                for i in m[-1][-1]:

                    if i == 'v':
                        mm += 'v'
                        m = [m[j] for j in range(len(m)-1)]
                        align_two_1(m, mm)

                    elif i == 'h':
                        mm += 'h'
                        m = [m[j][:-1] for j in range(len(m))]
                        align_two_1(m, mm)

        except IndexError:
            pass

    def generate_alignment(dna_0: str, dna_1: str, mm: str):

        d_0 = list(dna_0)
        d_1 = list(dna_1)

        res_0 = ''
        res_1 = ''
        indicator = ''
        n = 0

        # print(mm)
        for letter in mm:
            # print(letter, len(d_0), len(d_1), d_0, d_1)
            if letter == 'd':
                res_0 = d_0.pop() + res_0
                res_1 = d_1.pop() + res_1
            elif letter == 'h':
                res_0 = '-' + res_0
                res_1 = d_1.pop() + res_1
            elif letter == 'v':
                res_0 = d_0.pop() + res_0
                res_1 = '-' + res_1

        res_0 = ''.join(d_0) + res_0
        res_1 = ''.join(d_1) + res_1

        diff = len(d_0) - len(d_1)

        if diff > 0:
            res_1 = '-'*diff + res_1
        else:
            res_0 = '-'*(-diff) + res_0

        for i in range(min(len(res_1), len(res_0))):
            if res_0[i] == res_1[i]:
                indicator += '|'
                n += 1
            else:
                indicator += ' '

        all_res = res_0+'\n'+indicator+'\n'+res_1+'\n'

        # print(all_res)

        return n, all_res

    # f = open('score.csv', 'w')
    # for i in grid:
    #     for j in i:
    #         f.write(str(j))
    #         f.write(',')
    #     f.write('\n')
    # f.write('\n')
    # f.close()

    for i in range(20):
        for j in range(20):
            move[i][j] = ' '*(3-len(move[i][j]))+move[i][j]
    for i in move:
        print(i)
        print()
    align_two_1(move)
    print(set(way_of_alignment))
    sss = []
    for i in set(way_of_alignment):
        try:
            aaa = generate_alignment(dna_0, dna_1, i)
            sss.append(aaa)
        except:
            pass

    s = sorted(sss, key=lambda x: x[0], reverse=True)

    for i in s[:10]:
        print(i[1])

    return 0






def test():
    # a = DNA()
    #
    # a.set_name('name_A')
    # a.set_file_path('path')
    # a.set_seq('AGCTTAGCTAGCTAGCTGATCGATGCTAGCTG')
    #
    # print(a)
    #
    # a.rc()
    #
    # print(a, len(a))

    # linear_compare('AGCTTAGCTAAGCTTAGCTA', 'AGCTTCAGCTAGCTTAGCTA', [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10])

    # needleman('AGCTTAGCTAAGCTTAGCTA', 'AGCTTCAGACTACTTAGCTA',
    #           [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10], [-2]*20, [-2]*20)

    needleman('AATTAAAATTTTATTGACTT', 'CTAAATACTTTAACCAATAT',
              [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10], [-2]*20, [-2]*20)





def main():
    pass

if __name__ == '__main__':
    test()
