import pandas as pd
import os, pickle
import csv

with open('../data/pickle/pairwise_presence.unexpanded.pickle', 'rb') as f:
    df = pickle.load(f)

cells = df.columns.values

def makeCryptic(bounds, outdir):
    pairs5C = []
    for line in bounds:
        line = line.split('-')
        bound1 = line[0].split('_')
        bound2 = line[1].split('_')
        cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(bound1[0], bound1[1], bound1[2])
        cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(bound2[0], bound2[1], bound2[2])
        pairs5C.append([cb1, cb2])
    with open(outdir, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerows(pairs5C)
    return


for cell in cells:
    #need pos and neg pairfiles
    othercls = [k for k in cells if k != cell]
    otherBoundaries = df[othercls[0]] | df[othercls[1]] | df[othercls[2]] | df[othercls[3]]
    outdir = '../data/main/unexpanded/{}/'.format(cell)
    for ocl in othercls:
        pos = df[df[ocl] == True].index
        filename = '{}_posTADs.unexpanded.cryptic'.format(ocl)
        path = os.path.join(outdir, filename)
        makeCryptic(pos, path)
        remainingCls = [k for k in othercls if k != ocl]
        z = remainingCls    #alias
        # negative boundaries are those that do not exist in ocl
        # but do exist in at least one cell type from remainingCls
        neg = df[~df[ocl] & (df[z[0]] | df[z[1]] | df[z[2]])].index
        filename = '{}_negTADs.unexpanded.cryptic'.format(ocl)
        path = os.path.join(outdir, filename)
        makeCryptic(neg, path)


