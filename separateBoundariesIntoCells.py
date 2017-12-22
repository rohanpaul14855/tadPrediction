"""
Script to separate output files of getGlobalBoundaries
to separate files for each cell line and make them cryptic
"""

import csv
import os

files = [os.path.join('../data/single', k) for k in os.listdir('../data/single')
         if k.startswith('all')]
cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
path = '../data/single'
for file in files:
    pn = file.split('.')[2][-3:]        #positive or negative in filename
    pairs = list(csv.reader(open(file, 'r'), delimiter='\t'))
    for cl in cells:
        outname = '{}/{}_{}5CPairs.cryptic'.format(path, cl, pn)
        df = []
        for line in pairs:
            if line[0] == cl:
                df.append(line)
        out = []
        for line in df:
            cryptic_boundary1 = '5C_000_ENm000_FOR_000|hg19|chr{}:{}-{}'.format(line[1],
                                                                                int(float(line[2])),
                                                                                int(float(line[3])))
            cryptic_boundary2 = '5C_000_ENm000_REV_000|hg19|chr{}:{}-{}'.format(line[1],
                                                                                int(float(line[4])),
                                                                                int(float(line[5])))
            out.append([cryptic_boundary1, cryptic_boundary2])
            with open(outname, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(out)
