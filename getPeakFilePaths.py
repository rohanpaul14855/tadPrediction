import os

marks = '''H3K27ac
H3K27me3
H3K36me3
H3K4me1
H3K4me3
H3K9me3'''.split()

cells = [k for k in os.listdir('../data/21/peaks') if '.' not in k]
outdir = '../data/21/peaks/'

for cell in cells:
    cell_paths = {}
    paths = os.listdir('../data/21/peaks/{}'.format(cell))

    for mark in marks:
        for path in paths:
            if mark in path:
                cell_paths[mark] = os.path.join('/home/rpaul/looppred/data/21/peaks/{}/'.format(cell), path)


    outpath = os.path.join(outdir, '{}_paths.txt'.format(cell))
    with open(outpath, 'w') as f:
        for mark in sorted(cell_paths.keys()):
            f.write('{}\t{}'.format(mark, cell_paths[mark]))
            f.write('\n')

