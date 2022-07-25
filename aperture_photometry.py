#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.stats import sigma_clipped_stats
from matplotlib.patches import Circle
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.utils import calc_total_error
import argparse
import logging
import sys
from astropy.wcs import WCS

# In[2]:
def get_logger(namespace, level='DEBUG', logfile=None):
    logger = logging.getLogger(namespace)
    formatter = logging.Formatter('%(name)s [l %(lineno)d] - %(levelname)s - %(message)s')

    if logfile is None:
        handler = logging.StreamHandler(sys.stdout)

    else:
        handler = logging.FileHandler(logfile)

    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def make_cutouts(
        image_path: str,
        position,
        half_size
):
    data = fits.getdata(image_path)
    y_image_size, x_image_size = np.shape(data)
    x, y = position
    # logger.debug(f'{x},{y},{np.shape(data)}')
    if np.logical_and(x < half_size,y < half_size):
        cutout = data[0:y+half_size+1, 0:x+half_size+1]
        n_xpix = half_size-y
        n_ypix = half_size-x
        cutout = np.pad(cutout, ((n_ypix, 0), (n_xpix, 0)), 'constant')

    elif np.logical_and(x+half_size+1 > x_image_size, y+half_size+1 > y_image_size):
        cutout = data[y - half_size: y_image_size, x-half_size, x_image_size]
        n_xpix = (half_size+x+1) - x_image_size
        n_ypix = (half_size+y+1) - y_image_size
        cutout = np.pad(cutout, ((0, n_ypix), (0, n_xpix)), 'constant')

    elif y < half_size:
        logger.info(f'Cutout parameters are {y + half_size + 1}, {x - half_size}, {x + half_size + 1},{y_image_size},'
                    f'{x_image_size}')
        cutout = data[0:y + half_size + 1, x - half_size:x + half_size + 1]
        n_pix = half_size - y
        cutout = np.pad(cutout, ((n_pix, 0), (0, 0)), 'constant')

    elif y + half_size + 1 > y_image_size:
        cutout = data[y - half_size: y_image_size, x - half_size: x + half_size + 1]
        n_pix = (half_size + y + 1) - y_image_size
        cutout = np.pad(cutout, ((0, n_pix), (0, 0)), 'constant')

    elif x < half_size:
        cutout = data[y - half_size: y + half_size + 1, 0:x + half_size + 1]
        n_pix = half_size - x
        cutout = np.pad(cutout, ((0, 0), (n_pix, 0)), 'constant')
    elif x + half_size > x_image_size:
        cutout = data[y - half_size:y + half_size + 1, x - half_size:x_image_size]
        n_pix = (half_size + x + 1) - x_image_size
        cutout = np.pad(cutout, ((0, 0), (0, n_pix)), 'constant')
    else:
        cutout = data[y - half_size:y + half_size + 1, x - half_size:x + half_size + 1]
    return cutout


# In[73]:


def get_aperture_counts(diff_cutout, aper_diameter, bkg_in_diameter, bkg_out_diameter, x_offset = None, 
                    y_offset = None, gain=None, plot=False):

    #     w = WCS(header)
    #     x,y = w.all_world2pix(ra,dec,0)
    x, y = int(diff_cutout.shape[0] / 2), int(diff_cutout.shape[1] / 2)
    if x_offset is not None:
        x += x_offset
        y += y_offset
    if plot:
        fig, ax = plt.subplots()
        m, s = np.nanmean(diff_cutout), np.nanstd(diff_cutout)
        im = ax.imshow(diff_cutout, interpolation='nearest', cmap='gray',
                       vmin=m - s, vmax=m + 10 * s, origin='lower')
        # c = Circle(xy=(x_img, y_img),radius=15)

        c = Circle(xy=(x, y), radius=aper_diameter / 2)
        c1 = Circle(xy=(x, y), radius=bkg_in_diameter / 2)
        c2 = Circle(xy=(x, y), radius=bkg_out_diameter / 2)
        c.set_facecolor('none')
        c.set_edgecolor('red')
        c1.set_facecolor('none')
        c1.set_edgecolor('red')
        c2.set_facecolor('none')
        c2.set_edgecolor('red')
        ax.add_artist(c)
        ax.add_artist(c1)
        ax.add_artist(c2)
        ax.set_xlim(x - 30, x + 30)
        ax.set_ylim(y - 30, y + 30)
        plt.savefig(r'aper_phot.pdf',bbox_inches='tight')

    aperture = CircularAperture((x, y), r=aper_diameter)
    annulus_aperture = CircularAnnulus((x, y), r_in=bkg_in_diameter / 2, r_out=bkg_out_diameter / 2)

    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks.multiply(diff_cutout)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    bkg_mean, bkg_median, bkg_std = sigma_clipped_stats(annulus_data_1d, sigma=2)
    # print(bkg_mean, bkg_median)
    bkg = np.zeros(diff_cutout.shape) + bkg_median
    bkg_error = np.zeros(diff_cutout.shape) + bkg_std

    aperture_mask = aperture.to_mask(method='center')
    
    if gain is None:
        gain = 1
        print('Gain not provided, setting gain to 1, uncertainties will be incorrect (underestimated)')
        
    effective_gain = gain
    error = calc_total_error(diff_cutout, bkg_error, effective_gain)
    phot_table = aperture_photometry(diff_cutout - bkg, aperture, error=error)
    counts = phot_table['aperture_sum'][0]
    counts_err = phot_table['aperture_sum_err'][0]
    
    return counts, counts_err


def aper_photometry(imgname, x, y, zp, aper_diameter, bkg_in_diameter, bkg_out_diameter, gain=None, cutout_size=None, plot=False):
    x_int, y_int = int(x), int(y)
    x_offset, y_offset = x - x_int, y - y_int
    position=(x_int, y_int)
    
    if cutout_size is None:
        cutout_size = bkg_out+20
    
    cutout = make_cutouts(imgname,position,half_size=int(cutout_size))

    counts, countserr = get_aperture_counts(diff_cutout=cutout, aper_diameter=aper_diameter, bkg_in_diameter=bkg_in, bkg_out_diameter=bkg_out, x_offset=x_offset, y_offset=y_offset, gain=gain, plot=plot)

    mag = -2.5*np.log10(counts) + zp
    magunc = 1.086*countserr/counts

    return mag, magunc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("img",type=str,help="Image name")
    parser.add_argument("--ra",type=float,default=None,help="RA")
    parser.add_argument("--dec",type=float,default=None,help="Dec")
    parser.add_argument("--x",type=float,default=None,help="x coordinate")
    parser.add_argument("--y",type=float,default=None,help="y coordinate")
    parser.add_argument("--zp_key",type=str,default='TMC_ZP',help="zeropoint key in header")
    parser.add_argument("--gain_key",type=str,default='GAIN',help="gain key in header")
    parser.add_argument("--zp",type=float,default=None,help="zeropoint")
    parser.add_argument("--gain",type=float,default=None,help="gain")
    parser.add_argument("--aper",type=float,default=5,help="aperture diameter")
    parser.add_argument("--bkg_in",type=float,default=20,help="background inner annulus diameter")
    parser.add_argument("--bkg_out",type=float,default=30,help="background outer annulus diameter")
    parser.add_argument("--plot",action='store_true', help="Plot thumbnail with apertures")
    parser.add_argument("--cutout_size",type=float,default=None,help="Cutout size for display")


    args = parser.parse_args()

    logger = get_logger(__name__)
    

    imgname = args.img
    zp_key = args.zp_key
    gain_key = args.gain_key
    zp = args.zp
    cutout_size = args.cutout_size

    if zp is None:
        try:
            zp = float(fits.getval(imgname, zp_key))
        except KeyError:
            zp = 0
            logging.warning(f'Zeropoint not specified, or not found in header. Setting it to {zp}')

    x, y = args.x, args.y
    ra, dec = args.ra, args.dec
    if args.x is None:
        header = fits.getheader(imgname)
        wcs = WCS(header)
        if np.logical_or(ra is None, dec is None):
            err = 'Ra/Dec and x/y are not specified. Please specify either'
            logger.error(err)
            raise ValueError
        x, y = wcs.all_world2pix(ra,dec,0)

    logger.info(f'Setting x/y x : {x}, y:{y}')

    gain = args.gain

    if gain is None:
        try:
            gain = float(fits.getval(imgname, gain_key))

        except KeyError:
            gain = 1
            logger.warn(f'Gain not specified, or found in header. Setting to {gain}')


    aper_radius = args.aper
    bkg_in = args.bkg_in
    bkg_out = args.bkg_out
    mag, magunc = aper_photometry(imgname, x=x, y=y, zp=zp, aper_diameter=aper_radius, bkg_in_diameter=bkg_in, bkg_out_diameter=bkg_out, gain=gain, plot=args.plot, cutout_size=cutout_size)
    logger.info(f'Mag: {mag}, magerr: {magunc}')