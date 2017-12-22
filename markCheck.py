import os
import pandas as pd

cells = [k for k in os.listdir('../data/21/peaks/') if os.path.isdir('../data/21/peaks/' + k)]
marks = set()
for cell in cells:
    files = [k for k in os.listdir('../data/21/peaks/{}'.format(cell)) if k.endswith('.gz')]
    for file in files:
        marks.add(file.split('-')[1].split('.')[0])

master = {cell: {mark: False for mark in marks} for cell in cells}
for cell in cells:
    files = [k for k in os.listdir('../data/21/peaks/{}'.format(cell)) if k.endswith('.gz')]
    for file in files:
        mark = file.split('-')[1].split('.')[0]
        master[cell][mark] = True


master = pd.DataFrame(master)

master.to_csv('../data/21/peaks/mark_presence.csv', index=True)