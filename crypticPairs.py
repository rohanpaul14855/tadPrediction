#Converts pairfiles into format readable by sushmita paper code
#Outputs files with format 5C_###_ENm###_FOR_###|hg19|chr##:start-stop
#Seperates files by cell line

import csv
import os
import pandas as pd

def main(path=None):
    path = 'looplist/'
    if not os.path.isdir(path):
       os.makedirs(path)


    files = [k for k in os.listdir('looplist/') if k.endswith('TADs.txt') and not '5C' in k ]
    files = [os.path.join('looplist', k) for k in files]

    for file in files:
        cl = file.split('_')[0].split('/')[1]
        Pairs5C = []
        reader = csv.reader(open(file, 'r'), delimiter='\t')
        for line in reader:
            line = line[0].split('-')
            bound1 = line[0].split('_')
            bound2 = line[1].split('_')
            cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(bound1[0], int(float(bound1[1])), int(float(bound1[2])))
            cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(bound2[0], int(float(bound2[1])), int(float(bound2[2])))
            Pairs5C.append([cb1, cb2])
        outfile = 'looplist/{}_local_5CPairs.txt'.format(cl)
        writer = csv.writer(open(outfile, 'w'), delimiter='\t')
        print(len(Pairs5C))
        writer.writerows(Pairs5C)
main()
# positive_pairs = csv.reader(open('pos_pairwise.txt', 'r'), delimiter='\t')
# # positive_pairs = pd.DataFrame.from_csv('pos_pairwise.txt', sep='\t', index_col=None, header=None)
# for cl in ['HeLa', 'HMEC', 'HUVEC', 'K562', 'NHEK']:
#     outname = '{}/{}_positive5CPairs'.format(path, cl)
#     df = []
#     # df = list(positive_pairs[positive_pairs[positive_pairs.columns[0] == cl]])
#     for line in positive_pairs:
#         if line[0] == cl:
#             df.append(line)
#     out = []
#     for line in df:
#         cryptic_boundary1 = '5C_000_ENm000_FOR_000|hg19|chr{}:{}-{}'.format(line[1], int(float(line[2])), int(float(line[3])))
#         cryptic_boundary2 = '5C_000_ENm000_REV_000|hg19|chr{}:{}-{}'.format(line[1], int(float(line[4])), int(float(line[5])))
#         out.append([cryptic_boundary1, cryptic_boundary2])
#         # pd.DataFrame(out).to_csv(outname, sep='\t')
#         writer = csv.writer(open(outname,'w'), delimiter='\t')
#         writer.writerows(out)
#
# negative_pairs = csv.reader(open('neg_pairwise.txt', 'r'), delimiter='\t')
# # negative_pairs = pd.DataFrame.from_csv('neg_pairwise.txt', sep='\t', index_col=None, header=None)
# for cl in ['HeLa', 'HMEC', 'HUVEC', 'K562', 'NHEK']:
#     outname = '{}/{}_negative5CPairs'.format(path, cl)
#     df = []
#     for line in negative_pairs:
#         if line[0] == cl:
#             df.append(line)
#     # df = list(negative[negative_pairs.columns[0] == cl])
#     out = []
#     for line in df:
#         cryptic_boundary1 = '5C_000_ENm000_FOR_000|hg19|chr{}:{}-{}'.format(line[1], int(float(line[2])), int(float(line[3])))
#         cryptic_boundary2 = '5C_000_ENm000_REV_000|hg19|chr{}:{}-{}'.format(line[1], int(float(line[4])), int(float(line[5])))
#         out.append([cryptic_boundary1, cryptic_boundary2])
#         # pd.DataFrame(out).to_csv(outname, sep='\t')
#         writer = csv.writer(open(outname, 'w'), delimiter='\t')
#         writer.writerows(out)
