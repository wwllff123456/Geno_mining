import numpy as np



class DNA:

    def __init__(self):
        self.seq = 'N/A'
        self.name = 'N/A'
        self.file_path = 'N/A'

    def set_seq(self, seq: str):
        self.seq = seq

    def set_name(self, name: str):
        self.name = name

    def set_file_path(self, path: str):
        self.file_path = path

    def __str__(self):
        return 'DNA sequence {self.name} at {self.file_path}:\n{self.seq}'.format(self=self)



























def test():
    a = DNA()

    a.set_name('name_A')
    a.set_file_path('path')
    a.set_seq('AGCTTAGCTAGCTAGCTGATCGATGCTAGCTG')

    print(a)









def main():
    pass

if __name__ == '__main__':
    test()
