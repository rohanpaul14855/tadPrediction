import csv
import os

path = '../data/main/fixedNegatives'
cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']

for cell in cells:
    cpath = os.path.join(path, cell)
    fpath = os.path.join(cpath, '{}_fakeNegatives.boundaries'.format(cell))
    print(fpath)
    with open(fpath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        out = []
        for tad in reader:
            cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(tad[0], tad[1], tad[2])
            cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(tad[4], tad[5], tad[6])
            out.append([cb1, cb2])
        outfname = fpath.split('.')[2].split('/')[-1] +'.cryptic'
        outfname = os.path.join(path, cell, 'crypticBounaries', outfname)
        with open(outfname, 'w') as o:
            csv.writer(o, delimiter='\t').writerows(out)