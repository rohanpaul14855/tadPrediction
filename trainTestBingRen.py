"""
Creates training and testing cryptic pairs for each cell line/tissue
"""

import pickle
import os
import csv
import pandas as pd


def cryptify(pairs):
    """
    :param pairs: Pairs in the format chr_start_stop-chr_start_stop
    :return: 2 column list with cryptified pairs
    """
    out = []
    for tad in pairs:
        boundaries = tad.split('-')
        bound1 = boundaries[0].split('_')
        bound2 = boundaries[1].split('_')
        cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(bound1[0], bound1[1], bound1[2])
        cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(bound2[0], bound2[1], bound2[2])
        out.append([cb1, cb2])
    return out


with open('../data/21/pickle/pairwise_presence.unexpandedcont.pickle', 'rb') as f:
    presence = pickle.load(f)

traindir = '../data/21/cryptic/train'
testdir = '../data/21/cryptic/test'

cells = presence.columns.values

for cell in cells:

    # test
    pos_test_outfname = os.path.join(testdir, '{}_pos.cryptic'.format(cell))
    with open(pos_test_outfname, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        pos_tads = cryptify(presence[presence[cell]].index)     # true positive
        writer.writerows(pos_tads)

    neg_test_outfname = os.path.join(testdir, '{}_neg.cryptic'.format(cell))
    with open(neg_test_outfname, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        neg_tads = cryptify(presence[~presence[cell]].index)     # true negative
        writer.writerows(neg_tads)

    # train
    train_cell_df = presence.drop(cell, axis=1)
    train_out = os.path.join(traindir, cell)
    if not os.path.exists(train_out):
        os.makedirs(train_out)

    othercls = [k for k in cells if k != cell]
    for ocl in othercls:
        ocl_out = os.path.join(train_out, ocl)
        if not os.path.exists(ocl_out):
            os.makedirs(ocl_out)

        ocl_pos_outfname = os.path.join(ocl_out, '{}_pos.cryptic'.format(ocl))
        with open(ocl_pos_outfname, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            cryptic = cryptify(train_cell_df[ocl].index)
            writer.writerows(cryptic)

        ocl_neg_outfname = os.path.join(ocl_out, '{}_neg.cryptic'.format(ocl))
        with open(ocl_neg_outfname, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            remaining_cls = [k for k in othercls if k != ocl]
            keep_rows = []
            for _, row in train_cell_df.iterrows():
                keep_rows.append(True if row[ocl] and sum(row[remaining_cls]) != len(remaining_cls) else False)
            # True in ocl and False in at least one other cell ^^
            ocl_neg = train_cell_df[keep_rows].index
            cryptic = cryptify(ocl_neg)
            writer.writerows(cryptic)




    # neg train
    # pos test
    # neg test