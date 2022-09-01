import numpy as np
from astropy.io import fits
from glob import glob
import os

redo_list = ['ZTF22aazhflr','ZTF22abazrjk','ZTF22abbywgb','ZTF22abdbksc']
files = {}
for r in redo_list:
    files[r] = []
for r in ['bias','flat','arcs']:
    files[r] = []
for l in ls:
    h = fits.getheader(l)
    if h['object'] in files.keys():
        files[h['object']].append(l)
print(files)

for src in redo_list:
    os.mkdir(f'raw/{src}')
    for fname in files[src]:
        os.system(f'cp {fname} raw/{src}/')
    for caltype in ['bias','flat','arcs']:
        for fname in files[caltype]:
            os.system(f'cp {fname} raw/{src}/')