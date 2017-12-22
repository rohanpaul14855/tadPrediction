'''
Takes expanded TADs used to train the classifier and creates testing TADs from these
I.e. same TADs are used to test but with the held out cell line features
'''

import csv
import pickle
from os.path import join
import pandas as pd

with open('../data/pickle/pairwise_presence.unexpanded.pickle', 'rb') as f:
    pairs = pickle.load(f)

cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']

for cell in cells:
    ocl = [c for c in cells if c != cell]
    #Get positive pairs such that they exist in at least one other cell line
    #Only score prediction on pairs that can be predicted
    filtered_pairs = pairs[pairs[cell] & (pairs[ocl[0]] | pairs[ocl[1]] | pairs[ocl[2]] | pairs[ocl[3]])]
    pos_pairs = filtered_pairs[filtered_pairs[cell] == True].index
    neg_pairs = pairs[pairs[cell] == False].index

    def cryptify(tad):
        boundaries = tad.split('-')
        bound1 = boundaries[0].split('_')
        bound2 = boundaries[1].split('_')
        cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(bound1[0], bound1[1], bound1[2])
        cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(bound2[0], bound2[1], bound2[2])
        return [cb1, cb2]

    pos_pairs = list(map(cryptify, pos_pairs))
    neg_pairs = list(map(cryptify, neg_pairs))

    posoutfname = join('../data/main/unexpanded/test', '{}_pos_bobTestPairs.unexpanded.cryptic'.format(cell))
    negoutfname = join('../data/main/unexpanded/test', '{}_neg_bobTestPairs.unexpanded.cryptic'.format(cell))

    with open(posoutfname, 'w') as posout:
        csv.writer(posout, delimiter='\t').writerows(pos_pairs)

    with open(negoutfname, 'w') as negout:
        csv.writer(negout, delimiter='\t').writerows(neg_pairs)


