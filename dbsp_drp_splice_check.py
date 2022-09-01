from astropy.io import fits
from glob import glob
import numpy as np

ls = np.sort(ls)

for l in ls:
    header =fits.getheader(l)
    print('------------------')
    print(l.split('/')[-1].split('.fits')[0])
    print('------------------')
    print('Blue:',header['B_COADD'])
    print('Red:',header['R_COADD'])
