import time
from threading import Thread
from multiprocessing import Process


threads = 8

dna_path = 'example_data\\DNA.txt'

# output can be either txt or csv:
output_file = 'example_data\\output.csv'

# weight of match:
#         1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
weight = [1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 6, 8, 8, 8, 8]
# settings for Needleman-Wunsch:
match = 8
mis_match = -3
gap = -4

# guide is in forward sequence:
#           --------------------NGG
guide_20 = 'ACCATGCGAGTGTTGAAGTT'


# =================================================================


class Calculate:

    def __init__(self, dna_file, guide):
        self.t0 = time.time()

        self.dna_path = dna_file
        self.output = output_file
        self.guide20 = guide

        self.data = []        # list of dict, each dict repr a binging site
        self.p = [0]

    def extract(self):

        with open(self.dna_path, 'r') as d:

            dna = d.read().replace('\n', '')
            length = len(dna)

            dna_rc = dna.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]

            for i in range(21, len(dna) - 2):
                if dna[i:i+2] == 'GG':
                    temp_data = dict()
                    temp_data['from'] = i-21+1
                    temp_data['to'] = i+2+1
                    temp_data['seq'] = dna[i-21:i+2]
                    temp_data['rc'] = 'N'
                    self.data.append(temp_data)

                if dna_rc[i:i+2] == 'GG':
                    temp_data = dict()
                    temp_data['from'] = length - (i-21) + 1
                    temp_data['to'] = length - (i+2) + 1
                    temp_data['seq'] = dna[i-21:i+2]
                    temp_data['rc'] = 'Y'
                    self.data.append(temp_data)

    @staticmethod
    def needleman(dna_0: str, dna_1: str, m=match, mis=mis_match, g=gap) -> int:
        """Calculate the score_of_matching for two DNA sequences with Needleman-Wunsch algorithm"""

        grid = [[0 for _ in range(21)] for _ in range(21)]

        grid[0][1] = g
        grid[1][0] = g
        for i in range(2, 21):
            grid[0][i] = grid[0][i-1] + g
            grid[i][0] = grid[i-1][0] + g

        for i in range(20):
            for j in range(20):

                if dna_0[i] == dna_1[j]:
                    s = m
                else:
                    s = mis

                #               diagonal,           vertical,         horizontal:
                d, v, h = grid[i][j] + s, grid[i + 1][j] + g, grid[i][j + 1] + g

                grid[i+1][j+1] = max(d, v, h)

        return grid[-1][-1]

    def calc(self):

        class Progress(Thread):

            def __init__(self, now, total):
                super().__init__()
                self.now = now
                self.total = total

            def run(self):
                while True:
                    time.sleep(0.1)
                    print('{:.2%}%   {}/{}'.format(self.now[0]/self.total, self.now[0], self.total))
                    if self.now[0] >= self.total:
                        print('{:.2%}   {}/{}'.format(self.now[0] / self.total, self.now[0], self.total))
                        print('Done!')
                        break

        display = Progress(self.p, len(self.data)-1)
        display.start()

        # for perfect situation, the score would be:
        linear_p = 20
        weighted_p = sum(weight)
        needle_p = 20 * match

        for i, d in enumerate(self.data):

            self.p[0] = i

            seq = d['seq'][:20]

            linear = 0
            weighted = 0
            needle = self.needleman(seq, self.guide20)

            for j in range(20):
                if seq[j] == self.guide20[j]:
                    linear += 1
                    weighted += weight[j]

            d['linear'] = linear
            d['linear%'] = linear / linear_p

            d['weighted'] = weighted
            d['weighted%'] = weighted / weighted_p

            d['needleman'] = needle
            d['needleman%'] = needle / needle_p

        display.join()

    def sort_result(self, key='linear', reverse=False):
        """
        :param key: linear, werighted, needleman
        :param reverse: False: descending order, True: ascending order
        """
        self.data = sorted(self.data, key=lambda x: x[key], reverse=not reverse)

    def write_to_file(self):
        with open(self.output, 'w') as output:
            output.write('data for:,{}\n\n'
                         'from,to,rc,seq,linear,linear%,weighted,weighted%,'
                         'Needleman,Needleman%\n'.format(self.guide20))

            for each in self.data:
                output.write('{from},{to},{rc},{seq},'
                             '{linear},{linear%},{weighted},{weighted%},'
                             '{needleman},{needleman%}\n'.format(**each))


if __name__ == '__main__':
    t0 = time.time()

    a = Calculate(dna_path, guide_20)
    a.extract()

    a.calc()
    a.sort_result('weighted')
    a.write_to_file()

    t1 = time.time()

    print(round(t1-t0, 3))

