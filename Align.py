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

    def rc(self):
        logging.info('')
        self.seq = self.seq.upper()
        self.seq = self.seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
        self.seq = self.seq.upper()[::-1]
        logging.info('RC Done.')


    def __str__(self):
        return '\nDNA sequence {self.name} at {self.file_path}:\n{self.seq}\n'.format(self=self)



























def test():
    a = DNA()

    a.set_name('name_A')
    a.set_file_path('path')
    a.set_seq('AGCTTAGCTAGCTAGCTGATCGATGCTAGCTG')

    print(a)

    a.rc()

    print(a)









def main():
    pass

if __name__ == '__main__':
    test()
