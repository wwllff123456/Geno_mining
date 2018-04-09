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


def linear_compare(dna_0: str, dna_1: str, score_list: list, length=20) -> int:
    """Calculate the score_of_matching for two DNA sequences in a linear way"""

    if not len(dna_0) == len(dna_1) == len(score_list) == length:
        print('Not the same length in the two dna or the score list!!')
        logging.info('Length Error in the two dna or the score list: '
                     '\n[{}] {}\n[{}] {}\n[{}] {}\n'
                     'Length needs to be {}'
                     .format(len(dna_0), dna_0, len(dna_1), dna_1, len(score_list), score_list, length))
        logging.info('Exiting...')
        raise ValueError('Not the same length in the two dna or the score list!!')

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




    return 0


















def test():
    a = DNA()

    a.set_name('name_A')
    a.set_file_path('path')
    a.set_seq('AGCTTAGCTAGCTAGCTGATCGATGCTAGCTG')

    print(a)

    a.rc()

    print(a, len(a))

    linear_compare('AGCTTAGCTA', 'AGCTTCAGCT', [1,1,2,2,3,3,4,4,5,5])









def main():
    pass

if __name__ == '__main__':
    test()
