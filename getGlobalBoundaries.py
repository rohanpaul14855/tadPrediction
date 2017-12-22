# Generates positive and negative sample loops from .looplist.txt files
# provided by the rao paper

import csv
import os
import numpy as np
import pandas as pd
import random
import math
import pickle

cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
path = '/Users/rohanpaul/Dropbox/TadPrediction/data/looplist'
looplist = [k for k in os.listdir(path) if k.startswith('GSE')]


#reads chromosome sizes from file and returns size of specified chromosome
def getChromRange(chrom):
    chrom = 'chr{}'.format(chrom)
    sizes = csv.reader(open('../data/hg19.chrom.sizes', 'r'), delimiter='\t')
    for line in sizes:
        if line[0] == chrom:
            return int(line[1])
    else:
        exit('Invalid chromosome \"{}\", enter 1-22, X or Y'.format(chrom))


def genPosBoundaries():
    p = '../data/single/allpos.boundaries'
    outfile = csv.writer(open(p, 'w'), delimiter='\t')
    for loopfile in looplist:
        cl = loopfile.split('_')[1]
        pos_domains = []
        path = os.path.join('../data/looplist/', loopfile)
        infile = csv.reader(open(path, 'r'), delimiter='\t')
        next(infile)
        for line in infile:
            pos_domains.append([cl, line[0], line[1], line[2], line[4], line[5], 1])
        outfile.writerows(pos_domains)


def genNegBoundaries(sizes, path):
    num_boundaries = {'HeLa': 3094, 'HMEC': 5152, 'HUVEC': 3865, 'K562': 6057, 'NHEK': 4930}
    outfile = csv.writer(open('../data/single/allneg.boundaries', 'w'), delimiter='\t')

    for cell in cells:
        neg_domains = []
        loopfile = 'GSE63525_{}_HiCCUPS_looplist.txt'.format(cell)
        cl = loopfile.split('_')[1]
        abspath = os.path.join(path, loopfile)
        infile = csv.reader(open(abspath, 'r'), delimiter='\t')
        next(infile)
        for line in infile:
            rand = random.choice([0, 1])
            chrom, start_up, start_down, end_up, end_down = 0, 0, 0, 0, 0
            chrom = line[0]
            if rand == 0:     
                #pick both ends randomly

                start_possibilities = np.arange(5e4, getChromRange(chrom), 1e4)
                start_up = int(random.choice(start_possibilities))
                start_down = int(start_up + 1e4)
                while True:
                    chrom_idx = 0.5
                    if chrom.isdigit():
                        chrom_idx = int(chrom) - 1
                    elif chrom == 'X':
                        chrom_idx = 22
                    elif chrom == 'Y':
                        chrom_idx = 23
                    choices = np.arange(int(sizes[cl][chrom_idx][0]), int(sizes[cl][chrom_idx][1]), 1e4)
                    end_up = int(random.choice(choices))
                    # end_up = int(math.floor((end_up * std) + mean))
                    if end_up in range(int(line[4]) - 40000, int(line[4]) + 40000):
                        continue
                    else:
                        break
                end_down = int(end_up + 1e4)
            elif rand == 1:
                #fix one end
                chrom = line[0]
                start_up = int(line[1])
                start_down = int(line[2])
                ##no restrictions on size!?!
                choices = np.arange(start_down, int(sizes[cl][chrom_idx][1]), 1e4)
                if len(choices) == 0: continue
                end_up = int(random.choice(choices))
                ###########
                #TODO
                #Try removing the below while loop to see how prediction is affected
                #if there is no size restriction
                ##########
                while end_up > start_down + 2e6:
                    end_up = int(random.choice(choices))
                end_down = int(end_up + 1e4)
            neg_domains.append([cl, chrom, start_up, start_down, end_up, end_down, -1])
        outfile.writerows(neg_domains)



def size_dict(path):
    sizes = {}
    if not os.path.isfile('../data/pickle/sizes.pickle'):
        # print([k for k in os.listdir(path) if k.endswith('_looplist.txt')]); exit()
        for file in [k for k in os.listdir(path) if k.endswith('_looplist.txt')]:
            print(file)
            cl = file.split('_')[1]
            file = os.path.join(path, file)
            sizes[cl] = []
            for chrom in (list(range(1, 23)) + ['X', 'Y']):
                df = pd.DataFrame.from_csv(file, sep='\t', index_col=None)
                df = df[df['chr1'] == str(chrom)]
                sizes[cl].append((pd.DataFrame.min(df['x1']), pd.DataFrame.max(df['x1'])))
        with open('../data/pickle/sizes.pickle', 'wb') as f:
            pickle.dump(sizes, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open('../data/pickle/sizes.pickle', 'rb') as f:
            sizes = pickle.load(f)
    return sizes

# genPosBoundaries()
genNegBoundaries(size_dict('../data/looplist/'), '../data/looplist/')

