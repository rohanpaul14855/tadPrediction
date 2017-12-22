import os, csv
import pandas as pd
from pandas.core import index

import pickle

##################Collapses TADs and creates list of new expanded TADs

#####Increases size by 100kb on both sides of each boundary
#####TAD grows by 200kb
def expandTAD(TAD):
    tad_length = TAD[1][2] - TAD[0][1]
    bound1 = [TAD[0][0], TAD[0][1], TAD[0][2]]
    bound2 = [TAD[1][0], TAD[1][1], TAD[1][2]]
    return bound1, bound2


#####CHECK if tad1 contains tad2; true if yes, false otherwise
#####TAD2's boundaries must be completely inside the outer limits of TAD1
def contains(TAD1, TAD2):
    if TAD1[0][0] != TAD2[0][0]:    #check if on the same chromosome
        return False
    if TAD2[0][1] >= TAD1[0][1] and TAD2[1][2] <= TAD1[1][2]:
        return True
    return False


files = [k for k in os.listdir('../data/looplist') if k.endswith('looplist.txt') and not 'IMR90' in k]
files = [os.path.join('../data/looplist', k) for k in files]
cl = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']


master = [] #contains TADs from all cell lines
for file in files:
    cl = file.split('/')[-1].split('_')[1]
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        master += [[cl] + line[:6] + ['unvisited'] for line in reader] #exclude all the confidence measures



########################################################################################
# Does n^2 iterations through the master list                                          #
# For each observation, check all other boundaries to see if they fit in the first tad #
# If it does, marks it as present in that cell line                                    #
# Outputs a file with unique expanded TADs as rows and cell lines as columns           #
# Each entry in the table is True or False based on the presence in a cell line        #
########################################################################################

# if not os.path.isdir('../data/pickle/pairwise_presence.pickle'):
if True:
    master_presence, idx = [], []
    for i, line in enumerate(master):
        if line[-1] == 'visited':
            continue
        presence = {'HeLa': False, 'HMEC': False, 'NHEK': False, 'HUVEC': False, 'K562': False}
        cell = line[0]
        master[i][-1] = 'visited'
        tad1 = None
        try:
            tad1 = ([int(k) for k in line[1:4]], [int(k) for k in line[4:-1]])    #for chroms 1-22
        except:
            tad1 = ([line[1]] + [int(k) for k in line[2:4]], [line[4]] + [int(k) for k in line[5:-1]]) #for chroms X, Y
        tad1 = expandTAD(tad1)
        presence[cell] = True
        for j, check in enumerate(master):
            try:
                tad2 = ([int(k) for k in check[1:4]], [int(k) for k in check[4:-1]])
            except:
                tad2 = ([check[1]] + [int(k) for k in check[2:4]], [check[4]] + [int(k) for k in check[5:-1]])
            if contains(tad1, tad2):
                presence[check[0]] = True
                master[j][-1] = 'visited'

        pair = 'chr{0}_{1}_{2}-chr{0}_{3}_{4}'.format(line[1], line[2], line[3], line[5], line[6])
        master_presence.append(presence)
        idx.append(pair)
    master_presence = pd.DataFrame(master_presence, index=idx)
    with open('../data/pickle/pairwise_presence.unexpanded.pickle', 'wb') as f:
        pickle.dump(master_presence, f, pickle.HIGHEST_PROTOCOL)
else:
    with open('../data/pickle/pairwise_presence.pickle', 'rb') as f:
        df = pickle.load(f)

'''
with open('pickle/pairwise_presence(new).pickle', 'rb') as f:
    df = pickle.load(f)

cells = df.columns.values

for cell in cells:
    othercls = [k for k in cells if k != cell]
    bag = df[othercls[0]]| df[othercls[1]] | df[othercls[2]] | df[othercls[3]]
    print(cell, sum(bag))
    olap = sum(df[cell] & bag)
    totalBoundaries = sum(df[cell])
    # print('Overlapping Boundaries for {}: {}'.format(cell, olap))
    # print('Total Boundaries for {}: {}'.format(cell, totalBoundaries))
    # print('% overlap: {}'.format(olap*100/totalBoundaries))
    # print('\n')

'''

# print(cont)
one, two, three, four, five = 0, 0, 0, 0, 0
for i, row in master_presence.iterrows():
    numTrue = len([k for k in row if k == True])
    if numTrue == 1:
        one += 1
    if numTrue == 2:
        two += 1
    if numTrue == 3:
        three += 1
    if numTrue == 4:
        four += 1
    if numTrue == 5:
        five += 1

print('Total number of uniqe boundaries: {}'.format(master_presence.shape[0]))
print('Cell type specific: {}'.format(one))
print('Conserved across two cell types: {}'.format(two))
print('Conserved across three cell types: {}'.format(three))
print('Conserved across four cell types: {}'.format(four))
print('Conserved across five cell types: {}'.format(five))
