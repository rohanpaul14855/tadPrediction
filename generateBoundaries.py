# Generates positive and negative sample loops from .looplist.txt files
# provided by the rao et al. paper

import csv
import os
import numpy as np
import pandas as pd
import random
import math
import pickle
from bisect import bisect_left

cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
path = '../data/looplist/'
looplist = [k for k in os.listdir(path) if k.startswith('GSE')]

def getChromRange(chrom):
    """
    Get the length of the chromosome from hg19.chrom.sizes
    :param chrom: which chromosome
    :return: the length of the chromosome
    """
    chrom = 'chr{}'.format(chrom)
    sizes = csv.reader(open('../data/hg19.chrom.sizes', 'r'), delimiter='\t')
    for line in sizes:
        if line[0] == chrom:
            return int(line[1])
    else:
        exit('Invalid chromosome \"{}\", enter 1-22, X or Y'.format(chrom))

def getRandTad(cl, side, line):
    '''
    Generates locations for artificial tad boundaries 
    :param cl: The cell line in which the fake tad will exist
    :param side: 1: downstream boundary of the tad is fixed
    0: upstream boundary of the tad is fixed
    :param line: the line from .looplist file containing the coordinates of the real tad
    :return: the artificial tad [cl, chrom, start_up, start_down, end_up, end_down, -1]
    -1 for negative class
    '''
    tad_length = int(line[4]) - int(line[2])
    chrom = line[0]
    if side == 0:
        #side == 0, use upstream boundary of TAD
        end_up = int(line[1])   #fixed from real TAD
        end_down = int(line[2]) #fixed from real TAD
        start_down = end_up - tad_length    #fake boundary
        start_up = start_down - 10000       #fake boundary
        if start_up < 0 or start_down < 0:
            return None
        tosearch = [chrom, start_up, start_down]
        fixed = [chrom, end_up, end_down]
    else:
        #side == 1, use downstream boundary of TAD
        start_up = int(line[4])     #fixed from real TAD
        start_down = int(line[5])   #fixed from real TAD
        end_up = start_down + tad_length    #fake boundary
        end_down = end_up + 10000           #fake boundary
        tosearch = [chrom, end_up, end_down]
        fixed = [chrom, start_up, start_down]
    randtad = [str(k) for k in [cl, chrom, start_up, start_down, end_up, end_down]]
    if check(tosearch, fixed, side, cl):
        return randtad + [-1]
    else:
        return None

def check(tosearch, fixed, side, cl):
    '''
    Check whether the artificial boundary of this fake tad
    lies within another tad
    :param tad: The TAD whose boundaries are to be checked
    :return: True if boundary is not in another TAD, False otherwise
    '''
    chrom = tosearch[0]

    with open(os.path.join(path, 'GSE63525_{}_HiCCUPS_looplist.txt'.format(cl))) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            if line[0] == chrom:
                if int(line[1]) > int(tosearch[1]) + 1e5:
                    return True
                r1 = range(int(line[1]), int(line[5]) + 1)
                overlap = (int(tosearch[1]) in r1) or (int(tosearch[2]) in r1)
                match = False  #check if boundary is exact match
                if (fixed == line[:3]) or (fixed == line[3:6]):
                    match = True
                if overlap and not match:
                    return False

        else:
            #Iterated through the entire file and didn't find matches
            #Neither sides of the boundaries are in any of the TADs
            return True


def genPosBoundaries():
    outfile = csv.writer(open('pos_pairwise.txt', 'w'), delimiter='\t')
    for loopfile in looplist:
        cl = loopfile.split('_')[1]
        pos_domains = []
        path = os.path.join('looplist/', loopfile)
        infile = csv.reader(open(path, 'r'), delimiter='\t')
        next(infile)    #skip header
        for line in infile:
            pos_domains.append([cl, line[0], line[1], line[2], line[4], line[5], 1])
        outfile.writerows(pos_domains)


def genNegBoundaries(sizes, path):
    num_boundaries = {'HeLa': 3094, 'HMEC': 5152, 'HUVEC': 3865, 'K562': 6057, 'NHEK': 4930}
    for cell in cells:
        outfname = '../data/main/fixedNegatives/{}_fakeNegatives.boundaries'.format(cell)
        f = open(outfname, 'w')
        outfile = csv.writer(f, delimiter='\t')
        neg_domains = []
        loopfile = 'GSE63525_{}_HiCCUPS_looplist.txt'.format(cell)
        # cl = loopfile.split('_')[1]
        abspath = os.path.join(path, loopfile)
        infile = csv.reader(open(abspath, 'r'), delimiter='\t')
        next(infile)
        identifier = 0
        for line in infile:
            # rand = random.choice([0, 1])
            rand = 1    #only have tads fixed on one end
            chrom, start_up, start_down, end_up, end_down = 0, 0, 0, 0, 0
            chrom = line[0]
            if rand == 0:
                # pick both ends randomly
                start_possibilities = np.arange(5e4, getChromRange(chrom), 1e4)
                start_up = random.choice(start_possibilities)
                start_down = start_up + 1e4
                while True:
                    chrom_idx = 10000
                    if chrom.isdigit():
                        chrom_idx = int(chrom) - 1
                    elif chrom == 'X':
                        chrom_idx = 22
                    elif chrom == 'Y':
                        chrom_idx = 23
                    choices = np.arange(int(sizes[cl][chrom_idx][0]), int(sizes[cl][chrom_idx][1]), 1e4)
                    end_up = random.choice(choices)
                    # end_up = int(math.floor((end_up * std) + mean))
                    if end_up in range(int(line[4]) - 40000, int(line[4]) + 40000):
                        continue
                    else:
                        break
                end_down = end_up + 1e4
            else:
                # rand == 1; fix one end
                side = random.choice([0, 1])    #pick either side of the TAD
                randTad = getRandTad(cell, side, line)
                if randTad is None:             #check if random tad is valid
                    side = 1 if side == 0 else 0    #try with other side if so
                    randTad = getRandTad(cell, side, line)
                    if randTad is None: # check again
                        continue    # no negative for this tad
                # chrom = line[0]
                # start_up = int(line[1])
                # start_down = int(line[2])
                # choices = np.arange(int(sizes[cl][chrom_idx][0]), int(sizes[cl][chrom_idx][1]), 1e4)
                # end_up = random.choice(choices)
                # end_down = end_up + 1e4
            #write in the format [chr start1 stop1 chr start2 stop2]
            identifier += 1
            tad = ['chr' + randTad[1]] + randTad[2:4] + ['id{}'.format(identifier)] + ['chr' + randTad[1]] + randTad[4:-1] + ['id{}'.format(identifier)]
            print(tad)
            neg_domains.append(tad)
        outfile.writerows(neg_domains)
        f.close()




def size_dict(path):
    '''
    Uses sizes.pickle to get the range of each chromosome in which TADs exist
    This information is taken from the Rao et al paper
    :param path: path of loopfiles from rao et all
    :return: dictionary that contains the range for each chromosome in all cell lines
    '''
    sizes = {}
    if not os.path.isfile('../data/pickle/sizes.pickle'):
        for file in [k for k in os.listdir(path) if not k.startswith('.')]:
            cl = file.split('_')[1]
            file = os.path.join(path, file)
            sizes[cl] = []
            for chrom in (list(range(1, 23)) + ['X', 'Y']):
                df = pd.DataFrame.from_csv(file, sep='\t', index_col=None)
                df = df[df['chr1'] == str(chrom)]
                sizes[cl].append((pd.DataFrame.min(df['x1']), pd.DataFrame.max(df['x1'])))
        with open('pickle/sizes.pickle', 'wb') as f:
            pickle.dump(sizes, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open('../data/pickle/sizes.pickle', 'rb') as f:
            sizes = pickle.load(f)
    return sizes


# genPosBoundaries()
genNegBoundaries(size_dict('../data/looplist/'), '../data/looplist/')

#Call this bash script to turn it into a compatible file
# for i in *.boundaries; do name=$(echo $i | cut -d'.' -f1); cut -f1-4 $i > ${name}.left.bed; cut -f5-  $i > ${name}.right.bed; done
# cat a.sorted.bed | head | awk '{print $1"\t"$2"\t"$3"\t"$4}'
# mkdir Hela HMEC HUVEC NHEK K562
# for i in HeLa HMEC HUVEC NHEK K562; do mv ${i}_* $i; done

