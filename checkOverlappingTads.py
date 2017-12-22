import csv

def checkOverlap(interval1, interval2):
    c1, s1, e1 = interval1
    c2, s2, e2 = interval2

    if c1 != c2:
        return False

    if s1 < e2 and s2 < e1:
        #overlap
        return True
    else:
        return False

with open('../data/gm12878_bingren_tads.csv') as f:
    tads = list(csv.reader(f, delimiter=','))

for tad in tads:
    for tad2 in tads:
        if tad == tad2: continue
        if checkOverlap(tad, tad2):
            print(tad, tad2)
        exit()

