'''
Combines training features for each cell line

'''

import os
import csv
import numpy as np

cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
PATH = '../data/BWFeatures/main'

def getLocs(locs):
    header = np.array(['Locus'])
    for p in locs:
        if 'left' in p:
            with open(p, 'r') as f:
                tad_start = np.array(list(csv.reader(f, delimiter='\t')))[:, :3]
            break
    else:
        print('hwat')

    for p in locs:
        if 'right' in p:
            with open(p, 'r') as f:
                tad_end = np.array(list(csv.reader(f, delimiter='\t')))[:, :3]
            break
    else:
        print('hwat')
    tads = np.concatenate((tad_start, tad_end), axis=1).astype(object)
    print(tads.shape)
    func = lambda x: '-'.join(['_'.join(x[:3]), '_'.join(x[3:])])
    tads = np.array(list(map(func, tads)))
    print(tads.shape)
    tads = np.concatenate((header, tads), axis=0).reshape(-1, 1)
    return tads


for cell in cells:
    files = [os.path.join(PATH, cell, k) for k in os.listdir(os.path.join(PATH, cell))]
    for cell2 in [k for k in cells if k != cell]:

        pos_feats = [k for k in files if (cell2 in k) and ('pos' in k)]
        neg_feats = [k for k in files if (cell2 in k) and ('neg' in k)]
        # pos_locs = np.array(list(csv.reader(open(pos_feats[0], 'r'), delimiter='\t')))[:,:3]
        # neg_locs = np.array(list(csv.reader(open(neg_feats[0], 'r'), delimiter='\t')))[:,:3]
        pos_coordinates = getLocs(pos_feats)
        neg_coordinates = getLocs(neg_feats)

        for file in pos_feats:
            feature_name = file.split('.')[3]
            lr = file.split('.')[4]
            colname='{}_E'.format(feature_name) if lr == 'left' else '{}_P'.format(feature_name)
            with open(file, 'r') as f:
                feature_col = np.array(list(csv.reader(f, delimiter='\t')))[:,4]
                feature_col = np.concatenate(([[colname]], feature_col.reshape(-1, 1)), axis=0)

            #At each iteration, a feature column is concatenated to pos_coordinates
            pos_coordinates = np.concatenate((pos_coordinates, feature_col), axis=1)

        for file in neg_feats:
            feature_name = file.split('.')[3]
            lr = file.split('.')[4]
            colname='{}_E'.format(feature_name) if lr == 'left' else '{}_P'.format(feature_name)
            with open(file, 'r') as f:
                feature_col = np.array(list(csv.reader(f, delimiter='\t')))[:,4]
                feature_col = np.concatenate(([[colname]], feature_col.reshape(-1, 1)), axis=0)
            #At each iteration, a feature column is concatenated to neg_coordinates
            neg_coordinates = np.concatenate((neg_coordinates, feature_col), axis=1)
        #Now there is a pos coordinates and neg coordinates for a single cell line
        #Remove header for negative coordinates
        coordinates = np.concatenate((pos_coordinates, neg_coordinates[1:]), axis=0)
        outfname = os.path.join(PATH, cell, '{}.agg.features'.format(cell2))
        with open(outfname, 'w') as f:
            csv.writer(f, delimiter='\t').writerows(coordinates)
