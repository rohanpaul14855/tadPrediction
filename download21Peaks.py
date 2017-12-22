from bs4 import BeautifulSoup
import os
import requests
import urllib3


cell_codes = {'E003': 'hESC (H1)',
              'E005': 'Trophoblast-like Cell (TRO)',
              'E006': 'Mesenchymal Stem Cell (MSC)',
              'E007': 'Neural Progenitor Cell (NPC)',
              'E013': 'Mesendoderm (MES)',
              'E017': 'Fetal Lung Fibroblast (IMR90)',
              'E116': 'Lymphoblast (GM12878)',
              'E065': 'Aorta (AO)',
              'E066': 'Liver (LI)',
              'E071': 'Hippocampus (HC)',
              'E073': 'Dorsolateral Prefrontal Cortex (CO)',
              'E080': 'Adrenal (AD)',
              'E095': 'Left Ventricle (LV)',
              'E096': 'Lung (LG)',
              'E097': 'Ovary (OV)',
              'E098': 'Pancreas (PA)',
              'E100': 'Psoas (PO)',
              'E105': 'Right Ventricle (RV)',
              'E109': 'Small Bowel (SB)',
              'E113': 'Spleen (SX)',
              'None': 'Bladder (BL)'}



url = 'http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/'

html = requests.get(url).content
soup = BeautifulSoup(html, 'lxml')
for i, link in enumerate(soup.find_all('a')):
    fname = link.string.strip()
    if fname.endswith('.gz'):
        code = fname.split('-')[0]
        if code not in cell_codes.keys():
            continue
        mark = fname.split('-')[1].split('.')[0]
        download_url = url + fname
        cell_abbrv = cell_codes[code].split('(')[1][:-1]
        cell_name = cell_codes[code].split('(')[0][:-1]

        cell_peaks_path = '../data/21/peaks/{}/'.format(cell_abbrv)
        if not os.path.exists(cell_peaks_path):
            os.makedirs(cell_peaks_path)

        with open(cell_peaks_path + 'README.txt', 'w') as f:
            f.write(cell_name)

        outfile = cell_peaks_path + fname

        r = requests.get(download_url)
        with open(outfile, 'wb') as f:
            f.write(r.content)


